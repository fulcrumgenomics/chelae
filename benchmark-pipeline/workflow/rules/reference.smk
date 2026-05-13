# Reference FASTA handling. If `reference` in config is a URL, fetch it into
# results/reference/ (decompressing if .gz) and faidx it. Local paths are
# used as-is; the simulate rule always depends on `REFERENCE_LOCAL`.

rule fetch_reference:
    output:
        fa = REFERENCE_LOCAL,
        fai = REFERENCE_LOCAL + ".fai",
    params:
        url = config["reference"],
    resources:
        bench = 1  # blocked while a trim_one holds the bench pool
    log: "results/reference/fetch.log"
    shell:
        r"""
        mkdir -p "$(dirname {output.fa})"
        url="{params.url}"
        name=$(basename "$url")

        # Download the reference. Follow redirects; fail the rule on non-200.
        curl -sSfL "$url" -o "{output.fa}.dl" 2> {log}

        # Decompress if the URL ends in .gz; else move into place.
        if [[ "$name" == *.gz ]]; then
            gunzip -c "{output.fa}.dl" > "{output.fa}"
            rm "{output.fa}.dl"
        else
            mv "{output.fa}.dl" "{output.fa}"
        fi

        # Try to fetch the .fai alongside the source URL; fall back to faidx.
        if curl -sSfL "${{url}}.fai" -o "{output.fai}" 2>>{log}; then
            echo "fetched .fai sidecar" >> {log}
        else
            echo "no .fai sidecar; running samtools faidx" >> {log}
            samtools faidx "{output.fa}" 2>>{log}
        fi
        """
