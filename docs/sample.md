# Sample naming scheme (META-GP)

Metagenome samples (including TWIST samples from VIDRL) currently must follow a strict identifier naming scheme:

```
{SAMPLE_ID}__{NUCLEIC_ACID}__{SAMPLE_CATEGORY}__{REPEAT}_{ANYTHING_ELSE}
```

Note the double underscores (`__`) between the defined categories, single underscore (`_`) separate the defined components from any free text additional components (`ANYTHING_ELSE`) if necessary - these **must** be separated by single underscore from the rest of the identifier (`_`). 
The combined sample identifier must be unique across all samples that have ever been sequenced under this scheme.

where:

* `{SAMPLE_ID}` is the anonymised META-GP sample identifier, usually in the format: `DW-63-XXX` such as `DW-63-103` - this is the identifier of a single biological sample from a patient (or a single control like an environmental swab or sterile water control). If the same patient has a second biological sample taken - from the same site or a different site - this identifier should change. If a fresh water control is prepared this identifier should change. Do not use single underscores in this identifier, any components of the identifier should be separated by dashes (`-`), this component is **required**
* `{NUCLEIC_ACID}` is the library preparation nucleic acid and must be one of `DNA` or `RNA`, this component is **required**
* `{SAMPLE_CATEGORY}` must be one of the following categories`, this component is **required***:
  * `S` for a patient sample with unknown pathogen
  * `N` for a patient sample where we know they are negative for pathogens (for example we know they have cancer or autoimmune diagnosis)
  * `PS` for a patient sample with a known pathogen from orthogonal testing
  * `ENV` for environmental control swabs
  * `NTC` for negative template controls
  * `POS` for postitive control samples
* `{REPEAT}` should be `RPT1`, `RPT2`, ... where the number indicates the sequential repeat of a library constructed from the same biological sample (`SAMPLE_ID`), this component is **optional**.
* `{ANYTHING_ELSE}` is any other useful identifier text that may need to be doine


Important:

* Sample names must be anonymised and should **never** include linkable identifiers like Medipath numbers!
* There must **never** be spaces in the identifier!
* Do **not** include protected information like actual sample site or type like CSF or Blood, as this is protected patient information!
* If this naming scheme is not followed you will break all downstram applications with the potential for mis-interpretation of results and incorrect reporting to clinicians!

Examples:

A patient has two biological samples taken one from the brain (`DW-63-103`) and one from CSF (`DW-63-104`). DNA and RNA libraries are prepared for each biological sample. We do not know whether the patient has another diagnosis and we have no other orthogonal tests indicating a pathogen. This is not a repeat run of either biological sample. The sample is run on a large Illumina multiplex run where it is useful to put the index/sample from the input sample sheet into the name (1-4) just in case there are mixups. The sample names should therefore be:

```
DW-63-103__DNA__S_S1
DW-63-103__RNA__S_S2
DW-63-104__DNA__S_S3
DW-63-104__RNA__S_S4
```

Something happened and the libraries failed. In the next run a repeat of the brain and CSF samples are prepared. These are on different indices (18-21) in another large Illumina multiplex run. The sample names should therefore be:

```
DW-63-103__DNA__S__RPT1_S18
DW-63-103__RNA__S__RPT1_S19
DW-63-104__DNA__S__RPT1_S20
DW-63-104__RNA__S__RPT1_S21
```

A patient has a single biological sample taken from the brain (`DW-63-001`). RNA libraries are prepared for TWIST. We know this patient had an orthogonal testing positive for Ebola Virus (eesh). The sample names should therefore be:

```
DW-63-001__RNA__PS
```

Negative template controls from two distinct sterile water batches (`DW-63-WATER1`, `DW-63-WATER2`) and two environmental swabs (`DW-63-201E`, `DW-63-202E`) are prepared for a metagenomics run. DNA and RNA libraries are prepared for each. These are part of a large multiplex run on Illumina on indices 1-8. The sample names for these controls should therefore be:

```
DW-63-WATER1__DNA__NTC_S1
DW-63-WATER1__RNA__NTC_S2
DW-63-WATER2__DNA__NTC_S3
DW-63-WATER2__RNA__NTC_S4
DW-63-201E__DNA__ENV_S5
DW-63-201E__RNA__ENV_S6
DW-63-202E__DNA__ENV_S7
DW-63-202E__RNA__ENV_S8
```
