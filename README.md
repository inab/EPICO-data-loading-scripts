# NAME

BP::DCCLoader - EPICO data loading scripts

# SYNOPSIS

    git clone -b 20190711 https://github.com/inab/EPICO-data-model.git model
    carton exec epico-feature-loader.pl

# DESCRIPTION

BP::DCCLoader scripts are the data loading scripts created for
[BLUEPRINT project](https://blueprint-epigenome.eu).

        Then, you have to setup the project profile, either changing of copying L<epico-setup.ini|share/epico-setup.ini> file.

- Supported database paradigms (through \_**loaders**\_ configuration parameter in **mapper** section) are _relational_ (through \_**sql-dialect**\_ configuration
  parameter from **relational** section, subtypes _postgresql_, _mysql_ or _sqlite3_), _mongodb_ and _elasticsearch_ (these last ones have their own dedicated sections).
  On the dedicated section, uncomment and fill in only the parameters you need.
- You also have to setup several connection parameters, like \_**host**\_, \_**index-path**\_, \_**user**\_ and \_**pass**\_ on **dcc-loader** section in case you are using the BLUEPRINT data loading specific scripts.
- The usage of these scripts is:

    >     `perl epico-feature-loader.pl [-t] [-skipReactome] {iniFile} {cachindDir}`
    >
    >     `perl blueprint-dcc-loader.pl [-t|-tt] [-s|-ss] {iniFile} {cachingDir} [sdata|dnase|meth|pdna|rnaseq]`

    - `{inifile}` is the customized copy of [epico-setup.ini](https://metacpan.org/pod/share#epico-setup.ini) you have previously prepared.
    - `{cachingDir}` is a directory where the files to be processes are downloaded and cached. For primary data, those files are removed once they are processed.
    - Optional flags `-t` and `-tt` are used to put the scripts in test mode, so all the tasks are done but the insertion in the database.
    - Optional flags `-s` and `-ss` are to skip the insertion of primary data and its related specific analysis metadata.
    - Optional flag `-skipReactome` is usually used in combination with `-t` flag, in order to skip Reactome fetch and processing.

# METHODS

_(to be documented)_

# AUTHOR

José M. Fernández [https://github.com/jmfernandez](https://github.com/jmfernandez)

# COPYRIGHT

The library and programs were initially created several years ago for the
data management tasks in the
[BLUEPRINT project](http://www.blueprint-epigenome.eu/).

Copyright 2019- José M. Fernández & Barcelona Supercomputing Center (BSC)

# LICENSE

These programs and libraries are free software; you can redistribute them
and/or modify them under the Apache 2 terms.

# SEE ALSO
