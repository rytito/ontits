#!/usr/bin/env nextflow

/*
 * ontits -- Pipeline de Oxford Nanopore para amplicones ITS fúngicos
 *
 *   POD5 (sin basecalling) -> Dorado -> demux dual-index (cutadapt)
 *                          -> QC -> Emu (UNITE) -> phyloseq
 *
 * Derivado de ont16s. La diferencia central: la librería NO usa un kit de
 * barcoding de ONT sino primers dual-index estilo Illumina --
 * ITS1-F_KYO2 / ITS2_KYO2 con barcodes Golay de 12 nt (Toju et al. 2012,
 * esquema Earth Microbiome) -- así que `dorado demux` no aplica. El
 * demultiplexado usa cutadapt sobre la hoja assets/barcodes.tsv
 * (4 barcodes F x 4 barcodes R = 16 combinaciones).
 *
 * Escrito para la química Kit 14 (R10.4.1 / FLO-MIN114), librería ligada
 * sin barcodes nativos, con el basecalling DESACTIVADO en el secuenciador.
 *
 * Requiere Nextflow >= 26.04 (sintaxis estricta). La selección de etapa usa
 * --step, no -entry, que el parser estricto eliminó.
 */

/* ------------------------------------------------------------------ *
 *  Texto de ayuda
 * ------------------------------------------------------------------ */

def helpMessage() {
    log.info """
    ========================================================================
     ontits  v${workflow.manifest.version}
    ========================================================================

    Ejecución típica:

      nextflow run main.nf \\
          -profile slurm,conda \\
          -params-file conf/site.yml \\
          --step all \\
          --pod5_dir /ruta/a/la/corrida/pod5 \\
          --outdir results

    Etapas (--step, por defecto '${params.step}'):
      all         POD5 -> tablas de abundancia (+ phyloseq si se da --metadata)
      from_bam    --calls_bam <archivo>  omite el basecalling
      from_fastq  --fastq_dir <dir>      omite basecalling y demultiplexado

    Entradas:
      --pod5_dir      Directorio POD5. Con el basecalling desactivado MinKNOW
                      escribe 'pod5/', no 'pod5_pass/'.
      --calls_bam     BAM no alineado de un basecalling previo
      --fastq_dir     FASTQ por muestra (archivos llamados <SampleID>.fastq,
                      ya sin primers ni barcodes)
      --metadata      TSV con una columna SampleID; habilita la salida phyloseq
      --outdir        Directorio de salida [${params.outdir}]

    Basecalling:
      --model         Modelo de Dorado [${params.model}]
      --max_reads     Tope de lecturas (pruebas) [${params.max_reads ?: 'ninguno'}]

    Demultiplexado:
      --demux_mode         dual_index (Golay en primers, cutadapt) |
                           ont_kit (barcodes del kit + dorado demux)
                           [${params.demux_mode}]
      --barcodes           Hoja dual_index [assets/barcodes.tsv]
      --kit                Kit para ont_kit [${params.kit}]
      --ont_barcode_map    Hoja barcode->muestra para ont_kit
                           [assets/ont_barcode_map.tsv]
      --barcode_both_ends  Exigir barcode en ambos extremos (ont_kit) [${params.barcode_both_ends}]
      --fwd_primer         Primer forward [${params.fwd_primer}]
      --rev_primer         Primer reverse [${params.rev_primer}]
      --fwd_tag/--rev_tag  Tags universales ONT (ont_kit; ver docs/PRIMERS_ONT_RA.md)
      --demux_error_rate   Tasa de error de cutadapt [${params.demux_error_rate}]
      --demux_min_overlap  Solapamiento mínimo barcode+primer [${params.demux_min_overlap}]

    Filtrado y profundidad:
      --min_qscore    [${params.min_qscore}]
      --min_len       [${params.min_len}]
      --max_len       [${params.max_len}]
      --subsample     Lecturas por muestra, 0 = desactivado [${params.subsample}]
      --min_reads     Excluir muestras por debajo de esta profundidad [${params.min_reads}]

    Taxonomía:
      --emu_dbs       nombre:ruta[,nombre:ruta,...]
                      p. ej. 'unite:/db/emu_unite_fungi'

    Perfiles: slurm, vsc_wice, vsc_genius, conda, apptainer, test, debug

    Documentación completa: docs/USAGE.md (English: docs/USAGE.en.md)
    Lee docs/GOTCHAS.md antes de tu primera corrida real.
    """.stripIndent()
}

/* ------------------------------------------------------------------ *
 *  Funciones auxiliares
 * ------------------------------------------------------------------ */

// Los parámetros de la línea de comandos llegan como Strings; los de un
// config o params-file llegan tipados. Convierte antes de comparar números.
def asInt(v) {
    return v == null ? 0 : v.toString().toInteger()
}

// Complemento reverso con soporte IUPAC, para construir los adaptadores 3'
// de cutadapt a partir de los primers/tags declarados en params.
def revComp(s) {
    def M = ['A':'T','C':'G','G':'C','T':'A','U':'A','R':'Y','Y':'R','S':'S',
             'W':'W','K':'M','M':'K','B':'V','V':'B','D':'H','H':'D','N':'N']
    return s.toUpperCase().reverse().collect { M.get(it, 'N') }.join('')
}

// Dorado reconstruye el árbol de MinKNOW desde los metadatos del POD5:
//   <experiment>/<sample_id>/<run_id>/bam_pass/barcodeNN/*.bam
// sample_id suele ser literalmente "unknown". Nunca uses glob con ruta fija.
def barcodeOf(f) {
    def m = (f.toString() =~ /(barcode\d+|unclassified)/)
    return m ? m[0][1] : 'unknown'
}

// assets/ont_barcode_map.tsv -> [barcodeNN: SampleID]
def readBarcodeMap(f) {
    def m = [:]
    f.readLines().drop(1).each { line ->
        def bits = line.trim().split('\t')
        if( bits.size() >= 2 ) m[bits[0]] = bits[1]
    }
    return m
}

// "unite:/ruta,unite_all:/ruta" -> canal de [nombre, ruta]
def parseDatabases(spec) {
    if( !spec )
        error "--emu_dbs es obligatorio, p. ej. --emu_dbs 'unite:/db/emu_unite_fungi'"

    return Channel.from(spec.toString().split(','))
        .map { entry ->
            def bits = entry.trim().split(':')
            if( bits.size() != 2 )
                error "Entrada de --emu_dbs mal formada: '${entry}'. Se espera nombre:ruta"
            def dir = file(bits[1])
            if( !dir.exists() )
                error "Base de datos de Emu no encontrada: ${bits[1]}"
            tuple(bits[0], dir)
        }
}

/* ------------------------------------------------------------------ *
 *  Basecalling y demultiplexado
 * ------------------------------------------------------------------ */

process DORADO_BASECALL {
    label 'gpu'
    tag { pod5_dir.name }
    publishDir "${params.outdir}/01_basecall", mode: 'copy', pattern: '*.txt'

    input:
    path pod5_dir

    output:
    path 'calls.bam',    emit: bam
    path 'versions.txt', emit: versions

    script:
    def maxreads = params.max_reads ? "--max-reads ${params.max_reads}" : ''
    def bench    = params.batchsize_benchmarks ? "--batchsize-benchmarks-file ${params.batchsize_benchmarks}" : ''
    """
    nvidia-smi || { echo "ERROR: ninguna GPU visible en esta asignación" >&2; exit 1; }
    dorado --version > versions.txt 2>&1

    # --no-trim es esencial. El recorte automático de adaptadores de Dorado
    # puede comerse los extremos del amplicón, y los barcodes Golay + primers
    # que cutadapt necesita para demultiplexar están precisamente ahí.
    dorado basecaller ${params.model} ${pod5_dir} \\
        --no-trim \\
        --recursive \\
        --device cuda:all \\
        ${maxreads} ${bench} \\
        > calls.bam

    # Dorado puede salir con 0 tras una escritura truncada (la cuota de
    # disco es la causa habitual), así que verifica explícitamente. -u suprime
    # la advertencia "no targets in header", normal en BAM no alineado; solo
    # un bloque EOF ausente indica truncamiento real.
    samtools quickcheck -u -v calls.bam || {
        echo "ERROR: calls.bam está truncado -- revisa la cuota de disco" >&2
        exit 1
    }
    """
}

process BAM_TO_FASTQ {
    label 'cpu_low'

    input:
    path bam

    output:
    path 'all.raw.fastq', emit: fastq

    script:
    // -T '*' lleva las etiquetas de la corrida al campo de comentario del
    // FASTQ. El BAM se conserva como forma de archivo; esta es una
    // conversión de conveniencia para cutadapt.
    """
    samtools fastq -T '*' ${bam} > all.raw.fastq
    """
}

process CUTADAPT_DEMUX {
    label 'cpu_high'
    publishDir "${params.outdir}/02_demux", mode: 'copy',
        pattern: '{demux_stats.tsv,logs/*}'

    input:
    path fastq
    path barcodes

    output:
    path 'demux/*.fastq',    emit: fastqs
    path 'demux_stats.tsv',  emit: stats
    path 'logs/*',           emit: logs

    script:
    // Dos rondas de cutadapt: (1) orienta cada lectura (--revcomp) y la
    // asigna por barcodeF+primerF recortando el 5'; (2) asigna por
    // rc(barcodeR+primerR) recortando el 3'. Solo se aceptan lecturas con
    // AMBOS barcodes -- el equivalente de --barcode-both-ends: con el mismo
    // juego de 4 Golay en ambos extremos, el index hopping produciría
    // asignaciones cruzadas silenciosas.
    """
    demux_dual_barcodes.py \\
        --fastq ${fastq} \\
        --barcodes ${barcodes} \\
        --fwd-primer ${params.fwd_primer} \\
        --rev-primer ${params.rev_primer} \\
        --outdir demux \\
        --error-rate ${params.demux_error_rate} \\
        --min-overlap ${params.demux_min_overlap} \\
        --threads ${task.cpus}
    """
}

/* ------------------------------------------------------------------ *
 *  Demultiplexado modo ont_kit (SQK-RPB114.24 y afines)
 * ------------------------------------------------------------------ */

process DORADO_DEMUX {
    label 'cpu_high'
    publishDir "${params.outdir}/02_demux", mode: 'copy', pattern: '**summary*.txt'

    input:
    path bam

    output:
    path 'demux/**/*.bam',         emit: bams
    path 'demux/**summary*.txt',   emit: summary, optional: true

    script:
    def both = params.barcode_both_ends?.toString()?.toBoolean() ? '--barcode-both-ends' : ''
    """
    dorado demux ${bam} \\
        --output-dir demux \\
        --kit-name ${params.kit} \\
        ${both} \\
        --emit-summary \\
        --threads ${task.cpus}

    find demux -name '*.bam' | xargs -r samtools quickcheck -u -v
    echo "BAMs del demux verificados"
    """
}

process BAM_TO_FASTQ_BC {
    label 'cpu_low'
    tag { barcode }

    input:
    tuple val(barcode), path(bams)

    output:
    tuple val(barcode), path("${barcode}.raw.fastq"), emit: fastq

    script:
    """
    for f in ${bams}; do
        samtools fastq -T '*' \$f
    done > ${barcode}.raw.fastq
    """
}

process CUTADAPT_TRIM {
    label 'cpu_low'
    tag { sample }
    publishDir "${params.outdir}/02_demux/trim_logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.trim.fastq"), emit: fastq
    path "${sample}.cutadapt.log"

    script:
    // Tras dorado demux la lectura conserva [tagF][primerF]...rc([tagR][primerR]).
    // --revcomp orienta; --discard-untrimmed exige el primer F (las lecturas
    // sin sitio de primer no son amplicones).
    def a5 = "${params.fwd_tag}${params.fwd_primer}"
    def a3 = revComp("${params.rev_tag}${params.rev_primer}")
    """
    cutadapt -j ${task.cpus} \\
        -e ${params.demux_error_rate} -O 20 --revcomp \\
        --discard-untrimmed -m 50 \\
        -g ${a5} -a ${a3} \\
        -o ${sample}.trim.fastq ${fastq} > ${sample}.cutadapt.log
    """
}

/* ------------------------------------------------------------------ *
 *  QC de lecturas
 * ------------------------------------------------------------------ */

process CHOPPER_FILTER {
    label 'cpu_low'
    tag { sample }
    publishDir "${params.outdir}/03_filtered", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.fastq"), emit: fastq

    script:
    // A diferencia del 16S (ventana estrecha en ~1500 pb), la longitud de
    // ITS1 varía entre taxones fúngicos. El amplicón KYO2 de esta librería
    // mide ~300 pb (inserto ~260 pb ya sin primers ni barcodes); la ventana
    // deja margen a ambos lados para eliminar fragmentos y concatémeros sin
    // sesgar contra taxones de ITS más corto o más largo que el pico.
    """
    chopper -q ${params.min_qscore} \\
            --minlength ${params.min_len} \\
            --maxlength ${params.max_len} \\
            -i ${fastq} > ${sample}.fastq
    """
}

process SUBSAMPLE {
    label 'cpu_low'
    tag { sample }
    publishDir "${params.outdir}/04_subsampled", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.sub.fastq"), emit: fastq

    script:
    // Las abundancias relativas se estabilizan muy por debajo de 20k
    // lecturas, mientras que el tiempo de Emu escala con la profundidad.
    // La semilla fija lo hace reproducible.
    """
    seqkit sample -n ${params.subsample} -s ${params.seed} ${fastq} \\
        > ${sample}.sub.fastq
    """
}

process READ_STATS {
    label 'cpu_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path fastqs

    output:
    path 'read_stats.tsv'
    path 'reads_per_sample.tsv'

    script:
    """
    seqkit stats -a -T ${fastqs} > read_stats.tsv

    printf "sample\\treads\\n" > tmp.tsv
    for f in ${fastqs}; do
        s=\$(basename \$f | sed 's/\\..*//')
        n=\$(( \$(wc -l < \$f) / 4 ))
        printf "%s\\t%s\\n" "\$s" "\$n" >> tmp.tsv
    done
    ( head -1 tmp.tsv; tail -n +2 tmp.tsv | sort -k2,2nr ) > reads_per_sample.tsv
    rm tmp.tsv
    """
}

/* ------------------------------------------------------------------ *
 *  Perfilado taxonómico
 * ------------------------------------------------------------------ */

process EMU_ABUNDANCE {
    label 'cpu_med'
    tag { "${db_name}|${sample}" }
    publishDir { "${params.outdir}/05_emu/${db_name}" }, mode: 'copy'

    input:
    tuple val(sample), path(fastq), val(db_name), path(db_dir)

    output:
    tuple val(db_name), path("${sample}_rel-abundance.tsv"), emit: abundance

    script:
    // PYTHONNOUSERSITE evita que los paquetes de ~/.local eclipsen el entorno.
    """
    export PYTHONNOUSERSITE=1
    emu abundance ${fastq} \\
        --db ${db_dir} \\
        --type map-ont \\
        --keep-counts \\
        --threads ${task.cpus} \\
        --output-dir . \\
        --output-basename ${sample}
    """
}

process EMU_COMBINE {
    label 'cpu_low'
    tag { db_name }
    publishDir "${params.outdir}/06_tables", mode: 'copy', saveAs: { fn -> "${db_name}_${fn}" }

    input:
    tuple val(db_name), path(tsvs)

    output:
    tuple val(db_name), path('emu-combined-*.tsv'), emit: tables

    script:
    """
    export PYTHONNOUSERSITE=1
    mkdir -p combine && cp ${tsvs} combine/

    for rank in species genus family; do
        emu combine-outputs combine \$rank || true
        emu combine-outputs combine \$rank --counts || true
    done

    cp combine/emu-combined-*.tsv .
    """
}

/* ------------------------------------------------------------------ *
 *  phyloseq
 * ------------------------------------------------------------------ */

process BUILD_PHYLOSEQ {
    label 'cpu_low'
    tag { db_name }
    publishDir { "${params.outdir}/07_phyloseq/${db_name}" }, mode: 'copy'

    input:
    tuple val(db_name), path(tables)
    path metadata

    output:
    path 'phyloseq.rds',   emit: rds
    path '*.tsv',          emit: tsv,  optional: true
    path '*.pdf',          emit: pdf,  optional: true
    path 'sessionInfo.txt', emit: info, optional: true

    script:
    """
    COUNTS=\$(ls *combined-species-counts.tsv 2>/dev/null | head -1)
    if [ -z "\$COUNTS" ]; then
        echo "ERROR: no se encontró la tabla de conteos por especie" >&2
        exit 1
    fi
    build_phyloseq.R "\$COUNTS" ${metadata} . ${params.min_reads}
    """
}

/* ------------------------------------------------------------------ *
 *  Sub-workflows
 * ------------------------------------------------------------------ */

workflow DEMUX {
    take:
    ch_bam

    main:
    if( params.demux_mode == 'dual_index' ) {
        // Diseño original: barcodes Golay dual-index en los primers.
        BAM_TO_FASTQ(ch_bam)
        CUTADAPT_DEMUX(
            BAM_TO_FASTQ.out.fastq,
            file(params.barcodes, checkIfExists: true)
        )

        // demux/ contiene un FASTQ por SampleID de la hoja de barcodes. Las
        // lecturas sin ambos barcodes quedan en unassigned/ dentro de work/
        // y contabilizadas en demux_stats.tsv, pero fuera del análisis.
        ch_samples = CUTADAPT_DEMUX.out.fastqs
            .flatten()
            .map { f -> tuple(f.simpleName, f) }
    }
    else if( params.demux_mode == 'ont_kit' ) {
        // Diseño RA: primers KYO2 con tags universales + barcodes del kit
        // (SQK-RPB114.24). Demux por barcode nativo con dorado, luego
        // recorte de tags+primers con cutadapt. Ver docs/PRIMERS_ONT_RA.md.
        DORADO_DEMUX(ch_bam)

        ch_grouped = DORADO_DEMUX.out.bams
            .flatten()
            .map { f -> tuple(barcodeOf(f), f) }
            .filter { bc, f -> bc != 'unclassified' && bc != 'unknown' }
            .groupTuple()

        BAM_TO_FASTQ_BC(ch_grouped)

        // barcodeNN -> SampleID; los barcodes fuera de la hoja se descartan.
        bcmap = readBarcodeMap(file(params.ont_barcode_map, checkIfExists: true))
        ch_named = BAM_TO_FASTQ_BC.out.fastq
            .filter { bc, f -> bcmap.containsKey(bc) }
            .map { bc, f -> tuple(bcmap[bc], f) }

        CUTADAPT_TRIM(ch_named)
        ch_samples = CUTADAPT_TRIM.out.fastq
    }
    else {
        error "Valor de --demux_mode desconocido: '${params.demux_mode}'. Usa dual_index | ont_kit"
    }

    CHOPPER_FILTER(ch_samples)

    emit:
    fastq = CHOPPER_FILTER.out.fastq
}

workflow PROFILE {
    take:
    ch_fastq

    main:
    ch_in = asInt(params.subsample) > 0
          ? SUBSAMPLE(ch_fastq).fastq
          : ch_fastq

    READ_STATS(ch_in.map { s, f -> f }.collect())

    EMU_ABUNDANCE(ch_in.combine(parseDatabases(params.emu_dbs)))
    EMU_COMBINE(EMU_ABUNDANCE.out.abundance.groupTuple())

    if( params.metadata ) {
        BUILD_PHYLOSEQ(
            EMU_COMBINE.out.tables,
            file(params.metadata, checkIfExists: true)
        )
    }

    emit:
    tables = EMU_COMBINE.out.tables
}

/* ------------------------------------------------------------------ *
 *  Punto de entrada
 * ------------------------------------------------------------------ */

workflow {
    main:

    if( params.help ) {
        helpMessage()
    }
    else if( params.step == 'all' ) {
        if( !params.pod5_dir )
            error "--pod5_dir es obligatorio para --step all"
        ch_pod5 = Channel.fromPath(params.pod5_dir, type: 'dir', checkIfExists: true)
        DORADO_BASECALL(ch_pod5)
        DEMUX(DORADO_BASECALL.out.bam)
        PROFILE(DEMUX.out.fastq)
    }
    else if( params.step == 'from_bam' ) {
        if( !params.calls_bam )
            error "--calls_bam es obligatorio para --step from_bam"
        ch_bam = Channel.fromPath(params.calls_bam, checkIfExists: true)
        DEMUX(ch_bam)
        PROFILE(DEMUX.out.fastq)
    }
    else if( params.step == 'from_fastq' ) {
        if( !params.fastq_dir )
            error "--fastq_dir es obligatorio para --step from_fastq"
        ch_fastq = Channel
            .fromPath("${params.fastq_dir}/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: true)
            .map { f -> tuple(f.simpleName, f) }
        PROFILE(ch_fastq)
    }
    else {
        error "Valor de --step desconocido: '${params.step}'. Usa all | from_bam | from_fastq"
    }
}
