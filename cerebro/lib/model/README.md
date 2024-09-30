# `Cerebro` libraries and command-line interfaces

Main entry point to configure and deploy the `Cerebro` stack. Implements shared library components to keep each library and binary compilation as contained and small as possible.

* `cerebro`: command-line and library interface to user-oriented operations and deployment of the Cerebro stack

Imports and provides access to:

* `stack/watcher`: command-line and library interface to input file system upload watchers 
* `stack/workflow`: command-line and library interface to pipeline operations and outputs for Cerebro
* `stack/client`: command-line interface for a native client to the stack API

* `lib/model`: library interface to shared stack models between client and API

Does not import:

* `stack/server`: implementation and command-line runner of stack API imports `lib/model`