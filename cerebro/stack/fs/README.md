# `Cerebro` libraries and command-line interfaces

Main entry point to configure and deploy the `Cerebro` stack. Implements shared library components to keep each library and binary compilation as contained and small as possible.

* `cerebro`: command-line and library interface to user-oriented operations and deployment of the Cerebro stack

Imports and provides access to:

* `stack/fs`: command-line and library interface for handling SeaweedFS file storage/backup and integration Cerebro API
* `stack/report`: command-line and library interface for clinical reports from TOML or Cerebro API
* `stack/watcher`: command-line and library interface to input file system watchers and upload handlers
* `stack/workflow`: command-line and library interface to pipeline operations and outputs
* `stack/client`: command-line and library interface for a client to the Cerebro API
* `lib/model`: library interface to shared stack models between client and Cerebro API

Does not import:

* `stack/server`: implementation and command-line runner of stack API imports `lib/model`