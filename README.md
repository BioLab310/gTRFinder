# MTR_Detect

Release Date: March 20, 2025

Author

	~Changyong Yu (Northeastern University in CHINA)
	~Yuxing Qiao (Northeastern University in CHINA)
	~Yufan Chen (Northeastern University in CHINA)

1.Introduction

	MTR-Detect is a tool designed for accurate and efficient detection of tandem repeats (TRs) in genomic data. Based on the MiniDBG model, it utilizes De Bruijn Graph (DBG) and minimal edge extraction to improve detection in complex genomic regions. The tool performs well across different repeat lengths and sequencing errors, offering high precision and recall rates. MTR-Detect is faster and more reliable compared to existing tools, making it ideal for large-scale genomic analysis.

2.Test Data

	MTR-Detect employs the MinDBG structure to identify tandem repeats within the genome. The genomic sequence, provided as a .fa file, is utilized by MTR-Detect to generate in silico tandem repeats for testing purposes. The MTR-Detect tool can be invoked via the command-line interface to generate these simulated sequences.

3.Building Notes

Our project, built using the CMake build system, adheres to the following directory structure:

```
.
├── build 
├── CMakeLists.txt
├── config.h.in
├── include
├── lib
└── src
```

This codebase is developed in C++ and is compiled and executed on Linux-based operating systems. The prerequisites for building and running this software include the prior installation of CMake, GCC (GNU Compiler Collection), and Python version 3.9. The compilation process is as follows:

```b
mkdir build
cd build
cmake ..
make
```

We have uploaded a `build` directory, and we have successfully compiled the current files using the aforementioned procedure.

4.Usage Notes

1)Methods for TR identification in gene sequences

```
$ ./MTR_Detect -M -Sf <source file> [-Df <distination file>]

The `-M` command performs the identification of tandem repeats within the input gene sequence file (in `.fa` format). The identification results are stored in the `match_fasta/` directory, with the default naming convention inheriting the input filename. Alternatively, the user can specify the `-Df` option to rename the output files containing the identification results.
```

Running test

```
$ ./MTR_Detect -M -Sf 10_10.fasta -Df 10.fasta

The tool is invoked to analyze the tandem repeats within the `10_10.fasta` gene sequence file. The resulting tandem repeat sequences are saved in the `10.fasta` file.
```

2)Generation of \*in silico\* gene sequences

```
$ ./MTR_Detect -G [-L <len> -T <Repeat times> -Edr <edit distance ratios> -Sn <sequence_num -Cn <cell num> -OF <origin file> -Df <distinantion file]

The `-G` command is utilized to generate *in silico* gene sequences. Users have the option to input relevant parameters for the simulated sequences; if no parameters are provided, the generated sequence will be a random nucleotide sequence. The `-Df` option can be specified to designate the output filename; otherwise, the generated file will be automatically named `len_times`. The generated sequence file (in `.fasta` format) is saved in the `generate_fasta/` directory. Concurrently, a `.info` file containing the generation parameters is created in the same directory as the `.fasta` file.
```

Running test

```
$ ./MTR_Detect -G -L 10 -T 10 -edr 0.1 -Sn 10 -Cn 20 

The aforementioned command will generate a sequence file named `10_10.fasta`, containing ten gene sequences. Each gene sequence will feature 20 tandem repeats, with each repeat unit characterized by the following parameters: a repeat length of 10 nucleotides, a copy number of 10, an edit distance ratio of 0.1, and a randomly generated initial sequence. Concurrently, a `10_10.info` file will be created to store the detailed generation parameters of the sequences.
```


5.Parameter Settings

The format of a parameter in the MTR-Detect command line is a pair of strings, denoted here as `(-p, [q])`, `(-p, <q>)`, or `(-p)`. String `p` is the name of the parameter. String `q` is the value of the parameter input in the command line. `q` enclosed in square brackets `[]` represents an optional parameter. `q` enclosed in angle brackets `<>` represents a necessary parameter.

@parameter (-M)

	The -M parameter represents a method for invoking the tool to analyze sequences within a .fa formatted file.

@parameter (-Sf,<source file>)

	The -Sf parameter specifies the file containing the gene sequences to be analyzed. This file must be in .fa format.

@parameter (-Df,<distinnation file>)

	The -Df parameter specifies the filename for the output file containing the identified gene sequences. This file will be saved in the match_fasta/ directory.

@parameter (-G)

	The -G parameter represents a method for generating in silico tandem repeat sequence data. The generated sequence data is stored in the generate_fasta/ directory.

@parameter (-Of,<original file>)

	The -Of parameter is used for specifying the file containing the original gene sequence when generating in silico sequences.

@parameter (-T, <Repeat times/Range>)

	The T parameter specifies the copy number of the repeat unit within a generated tandem repeat sequence.
	
	The exact copy number can be set by providing a single integer value, such as -T 10.
	Alternatively, a range for the copy number can be defined by providing two integer values, such as -T 10 20.

@parameter (-L, <Repeat unit length/Range>)

	The L parameter specifies the length of the repeat unit within a generated tandem repeat sequence.
	
	The exact length can be set by providing a single integer value, such as -L 10.
	Alternatively, a range for the length can be defined by providing two integer values, such as -L 10 20.

@parameter (-Edr, <Edit distance ratio/Range>)

	The edr parameter specifies the proportion of the edit distance relative to a perfect tandem repeat sequence, with a valid range of [0, 1].
	
	The exact edit distance ratio can be set by providing a single floating-point value, such as -edr 0.1.
	Alternatively, a range for the edit distance ratio can be defined by providing two floating-point values, such as -edr 0 0.1.

@parameter (-Cn,<cell num>)

	The cn parameter specifies the number of tandem repeats to be inserted into a single gene sequence.

@parameter (-Sn,<sequence num>)

	The sn parameter specifies the number of gene sequences to be generated within the .fasta output file.

6.Format of Result

1)Identification TRs results

```
> sequence_num cell_num len interval_start interval_end others
repeated Sequences

Note: Due to thread concurrency, the sequence_number may not be strictly sequential.
```

2)Results of in silico sequence generation

```
.fasta
> sequence_num cell_num 

.info 
> sequence_num cell_num len interval_start interval_end repeated_times cell_sequence
```

7.License

	See LICENSE.txt

8.Contacts

	Please e-mail your feedback at cyyneu@126.com.



