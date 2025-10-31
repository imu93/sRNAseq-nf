// This may be over comented but its my first C code
// Just in case: https://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser
// Aim: collapse.c collapse identical FASTQ sequences keeping Rread Count (RC), and a mapping tsv
// build: gcc -O3 -march=native -pipe -DNDEBUG collapse.c -o collapse -lz
// usage: ./collapse -i reads.fastq.gz -f out.fastq.gz -M map.tsv

// Let's call the required libraries
#include <stdio.h>      // stdio for  printf/fprintf/fopen
#include <stdlib.h>     // C utils
#include <stdint.h>     // enteros de tamaño fijo
#include <string.h>     // permite manipular strings
#include <zlib.h>       // soporte gzip
#include <stdbool.h>    // logicos

// I will use the kseq and khash libraries to read fastq and call hashes
// The aim of kseq is to speed up reading hevy fastq, I do not whant to
// depend on my original siRmap module
#define KSEQ_INIT_STATIC
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)	// Kseq_init to read the gz fastq files, this is line by line

#include "khash.h"
// I will use hashes to build "dictionaries" (key=char*, val=uint64_t)
KHASH_MAP_INIT_STR(S, uint64_t)

// The hash will save seq quality and the counts per read instance
typedef struct {
    char    *seq;     // read sequence
    char    *qual;    // QUALITY of the firt read
    uint64_t read_count;     // Counts per identical read
} read_entry_t;

// Since in C I do not have list or vectors I will need to build one
typedef struct {
    read_entry_t *a;	// to point each new object to vec_t
    uint64_t n; // actual ellement
    uint64_t m; // total capacity of the vector
} read_vec_t;

// This functions will avoid the code die with no explanation
// static: just in this C
// void: no return
// die: die with msg
// and then just fprintf and exit
static void die(const char *msg){
    fprintf(stderr, "FATAL: %s\n", msg); // print message
    exit(1); // exit
}

// Now for memory with xmalloc
static void *xmalloc(size_t n) {
    void *p = malloc(n); // reserv mem
    if (!p)              // if there is not enough dead (Out of Memory)
        die("Out of memory");
    return p;            // if ok return p
}

// xalloc blocks and bites
static void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n, s); // start all the reserved memory at 0
    if (!p)
        die("Out of memory");
    return p;
}

// And to save new string copies in memory
static char *xstrdup(const char *s) {
    char *p = strdup(s); // duplicate string
    if (!p)
        die("Out of memory");
    return p;
}

static void readvec_reserve(read_vec_t *v, uint64_t cap){
    if (cap <= v->m) return;                          // ya tengo suficiente capacidad, nada que hacer
    uint64_t m = v->m ? v->m : 1024;                  // si está vacío, arranco con 1024 slots
    while (m < cap) m <<= 1;                          // duplico hasta cubrir 'cap'
    v->a = (read_entry_t*)realloc(v->a, m * sizeof(read_entry_t)); // agrando el buffer de entries
    if (!v->a) die("Out of memory (realloc)");        // por si falla el realloc, die
    v->m = m;                                         // actualizo la capacidad
}

static uint64_t readvec_push(read_vec_t *v, read_entry_t e){
    if (v->n == v->m) readvec_reserve(v, v->m ? v->m << 1 : 1024); // si no hay lugar, duplico o arranco en 1024
    v->a[v->n] = e;                                            // copio el entry en el siguiente slot
    return v->n++;                                             // retorno índice anterior y avanzo n
}

// The propose of this block is to set the hash for RID
// code  adapted from https://stackoverflow.com/questions/66764096/calculating-stdhash-using-different-compilers
static uint64_t fnv1a64(const char *s){
    const uint64_t FNV_PRIME = 1099511628211ULL;       // primo de FNV-1a
    uint64_t h = 1469598103934665603ULL;               // offset basis
    for (const unsigned char *p = (const unsigned char*)s; *p; ++p){
        h ^= (uint64_t)(*p);                           // xor con el byte actual
        h *= FNV_PRIME;                                // multiplicación por el primo
    }
    return h;                                          // hash 64-bit (luego lo recorto a 40 bits para el RID)
}

// function to set the out name
static bool ends_with_gz(const char *path){
    size_t L = strlen(path);                           // string length
    return (L >= 3 && strcmp(path + L - 3, ".gz") == 0); // sufix .gz??
}

// Now the main function
int main(int argc, char **argv){
    const char *input_fastq_path = NULL, *output_fastq_gz_path = NULL, *map_tsv_path = NULL;     // In- and Out- paths
    for (int i = 1; i < argc; ++i){                          // Now argsc like in optpasre or args
        if ((strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--in") == 0) && i + 1 < argc) input_fastq_path = argv[++i];
        else if ((strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--fastq") == 0) && i + 1 < argc) output_fastq_gz_path = argv[++i];
        else if ((strcmp(argv[i], "-M") == 0) && i + 1 < argc) map_tsv_path = argv[++i];
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0){
            fprintf(stderr, "usage: %s -i reads.fastq.gz -f out.fastq.gz -M map.tsv\n", argv[0]); // print help
            return 0;
        }
    }
    // Mandatory. -M is usefull for siRmap but for other utilities maybe not
    if (!input_fastq_path || !output_fastq_gz_path) die("need -i <fastq.gz> and -f <out.fastq.gz>");
    if (!ends_with_gz(input_fastq_path))  die("input FASTQ must be gzip (.gz)");
    if (!ends_with_gz(output_fastq_gz_path)) die("output FASTQ must be gzip (.gz)");

    gzFile in_gz = gzopen(input_fastq_path, "rb");
    if (!in_gz) die("cannot open input file");       // Open gz files

    // Now I will use Kseq to open the fastq
    kseq_t *seq = kseq_init(in_gz);
    if (!seq) die("kseq_init failed");

    // Open Khas
    // Let's build the dictionary seq id
    khash_t(S) *seq_index = kh_init(S);
    if (!seq_index) die("hash init failed");

    read_vec_t reads;
    memset(&reads, 0, sizeof(reads));                // beginner-style zero init (instead of compound literal)

    // conuter report of total and unique
    uint64_t total_reads = 0, unique_reads = 0;

    // pass 1: collapse by SEQ, keep first QUAL
    while (kseq_read(seq) >= 0){                    // Read block of 4 per read with kseq
        ++total_reads;                               // Save Total for my counter
        int absent = 0;                             // Internal flag
        // is my new sequence identic? if not this is a new instance
        khiter_t it = kh_get(S, seq_index, seq->seq.s);
        if (it == kh_end(seq_index)){
            // new unique
            read_entry_t e;
            memset(&e, 0, sizeof(e));
            // copy the new Seq
            e.seq = xstrdup(seq->seq.s);
            // And keep the first quality
            if (!seq->qual.s || seq->qual.l == 0)
                die("FASTQ record without quality string");
            if (seq->qual.l != seq->seq.l)
                die("sequence and quality length mismatch");
            e.qual = xstrdup(seq->qual.s);
            // I this the fisrt time this seq apears?
            e.read_count = 1;
            uint64_t idx = readvec_push(&reads, e);        // save in vec and take the index
            it = kh_put(S, seq_index, reads.a[idx].seq, &absent);
            kh_val(seq_index, it) = idx;
            ++unique_reads;
        } else {
            // seen before, retrieve the index and increase the RC of the read
            uint64_t idx = kh_val(seq_index, it);
            reads.a[idx].read_count += 1;
        }
    }

    // prepare outputs
    gzFile out_gz = gzopen(output_fastq_gz_path, "wb");   // salida gz obligatoria
    if (!out_gz) die("cannot open gz output");
    gzsetparams(out_gz, 4, Z_DEFAULT_STRATEGY);           // compression level 4 (fast)

    FILE *map_fp = NULL;                                  // archivo opcional para RID -> RC
    if (map_tsv_path){
        map_fp = fopen(map_tsv_path, "wb");               // TSV simple: rid \t rc (no gz)
        if (!map_fp) die("cannot open map output");
    }

    // pass 2: write FASTQ (first-seen order) and optional RID RC map
    for (uint64_t i = 0; i < reads.n; ++i){
        const read_entry_t *e = &reads.a[i];                 // tomo el entry i-ésimo (orden de primera aparición)
        uint64_t h = fnv1a64(e->seq) & 0xFFFFFFFFFFULL;      // low 40 bits                 // recorto a 10 hex
        char rid[64];                                        // buffer temporal para el ID
        snprintf(rid, sizeof(rid), "rid%010llx_%06llu",
                 (unsigned long long)h, (unsigned long long)(i + 1)); // id

        // Escribo registro FASTQ solo en gz (sin rama de texto plano)
        gzprintf(out_gz, "@%s\n%s\n+\n%s\n", rid, e->seq, e->qual);

        if (map_fp){
            fprintf(map_fp, "%s\t%llu\n", rid, (unsigned long long)e->read_count); // guardo el RID y su RC para downstream
        }
    }

    gzclose(out_gz);                         // cierro salida gz
    if (map_fp) fclose(map_fp);              // cierro TSV del mapa si existe

    // cleanup
    for (uint64_t i = 0; i < reads.n; ++i){              // libero cada string que dupliqué
        free(reads.a[i].seq);                            // free de la SEQ
        free(reads.a[i].qual);                           // free de la QUAL
    }
    // CAN NOT FORGET TO CLOSE ADN REMOVE HASH KSEQ AND INPUT
    free(reads.a);
    kh_destroy(S, seq_index);
    kseq_destroy(seq);
    gzclose(in_gz);
    // print small output
    fprintf(stderr, "done. total=%llu unique=%llu\n",
            (unsigned long long)total_reads, (unsigned long long)unique_reads);
    return 0;
}
