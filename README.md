# Kastor 
*(version 0.3.0)* 

Kastor is an error assessment program that uses aligned reference genome information to detect errors and suggest 
corrections in long read assemblies.

--- 

#### Modules
**Error detection**: Evaluate alignment information, n times, for potential errors above the error threshold. 

**Error adjustment**: Adjust detected putative errors based read data. Errors can only be adjusted 10 nucleotides 
upstream or downstream original position.  

**Correction**: Corrects the errors detected in the assembly by suggested corrections.
  
***

#### Requirements

**Packages**: Standard Python 3

--

**Files**: 

*Error detection*
- Ref mpileup file: Reference genomes aligned to draft assembly in a samtools-formated mpileup file.
- Optional: Supplementary ref mpileup file for low depth region. This file can only be declared with an --low 
flag

*Error Adjustment*
- Read mpileup file: Original reads aligned to draft assembly
- Optional: Additional mpileup  files. These can include pre-polished assemblies aligned to the draft assembly.

*Correction*
- Draft assembly
- Kastor error output file: This file needs to be formatted to post-adjustment format; substitutions cannot contain 
more one suggested output. *(Not yet implemented)* Set --clean to output pre-adjustment as post-adjustment format;
one caveat is that all substitution possibilities will not be considered if adjusted later with the adjustment 
module.  

---

##### Example Run

Full workflow
```bash
# Full pipeline
/path/to/kastor.py Ref.mpileup --reads Read.mpileup --draft Assm.fa [OPTIONS]

# Multiple adjustment files (Adjustment files are used sequentially)
/path/to/kastor.py Ref.mpileup --reads Read1.mpileup Read2.mpileup --draft Assm.fa [OPTIONS]

# Use a supplemental reference mpileup file
/path/to/kastor.py Ref.mpileup --reads Read.mpileup --draft Assm.fa --low supp_ref.mpileup [OPTIONS] 
```

Single/Double module runs 

<sub>*NOTE: Currently, files are read as string to be read and opened, thus string order still needs to be maintained* </sub>
```bash
# Error detection
/path/to/kastor.py Ref.mpileup --module Detect [OPTIONS]

# Error adjustment
/path/to/kastor.py Ref.mpileup --reads Read.mpileup --module Adjust --errors kastor_out.err [OPTIONS]

# Correction (Post-adjustment format ONLY)
/path/to/kastor.py Ref.mpileup --draft Draft.fa --module Correct --errors kastor_out.err [OPTIONS]
```
---

#### Options

|Abbrev|Option           |Requirements|Default|Description                                                            |
|------|-----------------|------------|-------|-----------------------------------------------------------------------|
|`-h`  | `--help`        |            |       | Shows this help message and exit                                      |
|**Thresholds Options**                                                                                               |
|`-d`  |`--depth`        | *int*      | 5     | Minimum depth required for detection                                  |
|`-i`  |`--indel`        | *float*    | 0.99  | Indel error threshold for detection                                   |
|`-s`  |`--sub`          | *float*    | 1.0   | Substitution error threshold for detection                            |
|      |`--low_threshold`| *float*    | 0.99  | Indel error threshold for supplementary reference alignment file      |
|`-a`  |`--adjust`       | *float*    | 0.25  | Adjustment error threshold                                            |
|**Others**                                                                                                           |
|`-l`  |`--low`          | *file*     |       | Enter a supplementary reference alignment file for low depth sequences|
|`-m`  |`--module`       | *str*      | "All" | Select module to be run                                               |
|`-p`  |`--passes`       | *int*      | 4     | Go through reference data n times, increasing search window each pass |
|`-o`  |`--out`          | *str*      | "out" | supply a prefix for output files                                      |
|      |`--errors`       | *file*     |       | Bypass error detection by entering an alternative error file          |
|**Misc**                                                                                                             |
|      |`--silent`       |            | False | Turn off outputs to stdout                                            |
|`-e`  |`--extraInfo`    |            | False | Print all module outputs and extra diagnostic data (Needs: draft)     | 
|`-v`  |`--version`      | *int*      |       | Software version number                                               |

---

#### Citation 
*To be announced*
 
#### License
    Kastor: reference-based error detection and correction software for
    long-read assemblies
    Copyright (C) 2024  Janet Lorv

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

