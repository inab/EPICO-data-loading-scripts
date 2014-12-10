BLUEPRINT-dcc-loading-scripts
=============================

In order to use BLUEPRINT DCC database loading scripts you have to clone
the repository:

	git clone --recurse-submodules https://github.com/inab/BLUEPRINT-dcc-loading-scripts.git
	
and the update to the latest changes is with:

	git pull --recurse-submodules
	git submodule update --recursive

Then, you have to setup the project profile, either changing of copying [blueprint-setup.ini](blueprint-setup.ini) file.

* Supported database paradigms (through _**loaders**_ configuration parameter in **mapper** section) are *relational* (through _**sql-dialect**_ configuration
  parameter from **relational** section, subtypes *postgresql*, *mysql* or *sqlite3*), *mongodb* and *elasticsearch* (these last ones have their own dedicated sections).
  On the dedicated section, uncomment and fill in only the parameters you need.

* You also have to setup the connection parameters _**host**_, _**index-path**_, _**user**_ and _**pass**_ to BLUEPRINT data repository on **dcc-loader** section.

