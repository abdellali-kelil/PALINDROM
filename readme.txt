PALINDROM Algorithm v1.0
________________________________________________________
Abdellali Kelil, PhD
Reseach Associate and Senior Bioinformatician
Sidhu Laboratory - Donnelly Centre - University of Toronto
160 College Street Room 810 - Toronto, Ontario CANADA M5S 3E1
abdellali.kelil@utoronto.ca

Description
________________________________________________________
PALINDROM Algorithm C/C++ implementation to find all palindromic sequences of length L with
defined nucleotide bases (no ambiguous bases) in a long DNA sequence (e.g. a chromosome).

Output
________________________________________________________

There are two output files
	1. DATA\\Times.dat
	2. DATA\\Palindromes.dat

The first output file has two columns:
	1. Column one is the alindromic sequences present in the long input sequence sorted alphabetically.
	2. Column two includes a comma-delimited string of the positions (one-based coordinate, sorted in
       ascending order) of the occurrences of the corresponding palindrome in the input sequence.

The second output file has three sections:
	1. Input sequences : gives input sequence information
	2. Find palindromes in segment : give information on output palindromes
	3. Execution time : gives average and stadard deviation of benchmark execution time

Execution
________________________________________________________
PALINDROM Algorithm was developped using VS C++ using standard libraries.
Open the project on VS C++ 2019 (free community edition), and execute as it is!!!
You can also, compile the project using VS C++, and run with the terminal.

Usage
________________________________________________________
 Usage : PALINDROM [file] [start] [length] [palindromes] [executions]
 Arguments :
 [file]        : Fasta input file (one or more sequences allowed)
 [start]       : segment start position (start the search at this position)
 [length]      : segment length (search within this kength)
 [palindromes] : palindromes length
 [executions]  : number of executions (evaluates execution time)
