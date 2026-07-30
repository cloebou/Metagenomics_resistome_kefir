process HOST_DB {

    tag "host_db"
    publishDir "data/host", mode: 'copy'

    conda "conda-forge::wget=1.21.4 conda-forge::unzip=6.0"

    input:
        val link_db

    output:
        path "host/*", type: 'dir'

    script:
    
    """
    mkdir -p host
    cd host

    for url in ${link_db.join(' ')}; do
        fname=\$(basename \$url)
        species=\${fname%.zip}

        wget -q \$url
        unzip -q \$fname
        rm -f \$fname
    done
    """
}
