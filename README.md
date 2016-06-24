EPICO-data-loading-scripts
=============================

In order to use both EPICO data loading scripts and BLUEPRINT ones you have to clone
the repository:

	git clone --recurse-submodules https://github.com/inab/EPICO-data-loading-scripts.git
	
and, in case you already have a copy, the update to the latest changes is with:

	git pull --recurse-submodules
	git submodule update --recursive

Then, you have to setup the project profile, either changing of copying [epico-setup.ini](epico-setup.ini) file.

* Supported database paradigms (through _**loaders**_ configuration parameter in **mapper** section) are *relational* (through _**sql-dialect**_ configuration
  parameter from **relational** section, subtypes *postgresql*, *mysql* or *sqlite3*), *mongodb* and *elasticsearch* (these last ones have their own dedicated sections).
  On the dedicated section, uncomment and fill in only the parameters you need.

* You also have to setup several connection parameters, like _**host**_, _**index-path**_, _**user**_ and _**pass**_ on **dcc-loader** section in case you are using the BLUEPRINT data loading specific scripts.

* The usage of these scripts is:

	```
	perl epico-feature-loader.pl [-t] [-skipReactome] {iniFile} {cachindDir}
	```
	
	```
	perl blueprint-dcc-loader.pl [-t|-tt] [-s|-ss] {iniFile} {cachingDir} [sdata|dnase|meth|pdna|rnaseq]
	```
	
  * `{inifile}` is the customized copy of [epico-setup.ini](epico-setup.ini) you have previously prepared.
  
  * `{cachingDir}` is a directory where the files to be processes are downloaded and cached. For primary data, those files are removed once they are processed.
  
  * Optional flags `-t` and `-tt` are used to put the scripts in test mode, so all the tasks are done but the insertion in the database.
  
  * Optional flags `-s` and `-ss` are to skip the insertion of primary data and its related specific analysis metadata.
  
  * Optional flag `-skipReactome` is usually used in combination with `-t` flag, in order to skip Reactome fetch and processing.