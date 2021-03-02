//---------------------------------------------------------------------------
//
#include "Utils.hpp"

//---------------------------------------------------------------------------
// Read fasta file format
__int64 ReadFastaFile(char *inp, __int64 &NBR, __int64 *&LEN, char **&SEQ)
{
	// Buffer
	char *ptr;

	//////////////////////////////////////////////
	/////////////////// STEP 1 ///////////////////
	//////////////////////////////////////////////

	// Verbose
	printf_s("STEP 1 : Fasta file size\n");

	// init lengths and size
	__int64 ln1 = 0; __int64 sze = 0;

	// Verbose
	__int64 nb0 = 0; __int64 nb1 = 0; __int64 nb2 = 0;

	// positions
	__int64 ps1 = 0; __int64 ps2 = 0;

	// init lengths
	NBR = 0; LEN = new __int64[1];

	// Get input filename
	FILE *file1 = OpenFile(inp, "rb");

	// get file size in bytes
	_fseeki64(file1, 0, SEEK_END); sze = _ftelli64(file1);

	// not end of file
	while (!feof(file1))
	{
		// back to last header
		_fseeki64(file1, ps1, SEEK_SET);

		// wrtite to file
		if (fread(BUFFER0, sizeof(char), MAXBUF, file1) == 0) break;

		// get buffer
		char *buf = BUFFER0;

		// read seq
		while (buf != NULL)
		{
			/////////////////// check new line ///////////////////

			// check new line
			if (*buf == '\n' || *buf == '\r')
			{
				// next
				buf += 1;

				// check read
				if (buf - BUFFER0 == MAXBUF) break;

				// check end of buffer
				if (ps1 + buf - BUFFER0 == sze) break;

				// next
				continue;
			}

			////////////////// read from buffer //////////////////

			// pointer
			ptr = buf;

			// length
			LEN1 = 0;

			// scroll
			while (true)
			{
				// update sequence
				BUFFER1[LEN1] = *ptr;

				// check end of sequence										// next
				if (*ptr == '\n' || *ptr == '\0' || LEN1+1 == MINBUF-1) break;	ptr++;

				// length
				LEN1 += 1;
			}

			// end sequence and next
			BUFFER1[LEN1+1] = '\0'; ptr += 1;

			//////////////// check end of buffer /////////////////

			// check read
			if (ptr - BUFFER0 > MAXBUF)
			{
				// correct length
				LEN1 = BUFFER0 - buf + MAXBUF;

				// trunk buffer
				BUFFER1[LEN1] = '\0';

				// update pointer
				ptr = buf + LEN1;
			}

			// check end of buffer
			if (ps1 + ptr - BUFFER0 > sze)
			{
				// correct length
				LEN1 = (sze - ps1) - (buf - BUFFER0);

				// trunk buffer
				BUFFER1[LEN1] = '\0';

				// update pointer
				ptr = buf + LEN1;
			}

			//////////////// read header/sequence ////////////////

			// Check if seq id
			if (BUFFER1[0] == '>')
			{
				// check end of line
				if (strchr(BUFFER1, '\n') == NULL) break;

				// Inc seq
				NBR += 1;

				// Init len
				LEN = Realloc(LEN, NBR-1, NBR); LEN[NBR-1] = 0;
			}
			else
			{
				// update lengths
				nb0 += LEN1; LEN[NBR-1] += LEN1;
			}

			// pointer
			buf = ptr;

			// check end of buffer
			if (buf - BUFFER0 >= MAXBUF) break;

			// check end of buffer
			if (ps1 + buf - BUFFER0 >= sze) break;
		}

		// update position
		ps1 += buf - BUFFER0;

		// Verbose
		printf_s("length : %15s\r", format(nb0, ','));
	}

	fclose(file1);

	// Verbose
	printf_s("\n");

	//////////////////////////////////////////////
	/////////////////// STEP 2 ///////////////////
	//////////////////////////////////////////////

	// Verbose
	printf_s("STEP 2 : Fasta file read\n");

	// Init sequences
	SEQ = new char*[NBR]; for (__int64 i=0; i<NBR; i++) SEQ[i] = new char[LEN[i]+1];

	// Get input filename
	FILE *file2 = OpenFile(inp, "rb");

	// scroll sequences
	while (!feof(file2))
	{
		// back to last header
		_fseeki64(file2, ps2, SEEK_SET);

		// wrtite to file
		if (fread(BUFFER0, sizeof(char), MAXBUF, file2) == 0) break;

		// get buffer
		char *buf = BUFFER0;

		// read seq
		while (buf != NULL)
		{
			/////////////////// check new line ///////////////////

			// check new line
			if (*buf == '\n' || *buf == '\r')
			{
				// next
				buf += 1;

				// check read
				if (buf - BUFFER0 == MAXBUF) break;

				// check end of buffer
				if (ps2 + buf - BUFFER0 == sze) break;

				// next
				continue;
			}

			////////////////// read from buffer //////////////////

			// pointer
			ptr = buf;

			// length
			LEN1 = 0;

			// scroll
			while (true)
			{
				// update sequence
				BUFFER1[LEN1] = *ptr;

				// check end of sequence										// next
				if (*ptr == '\n' || *ptr == '\0' || LEN1+1 == MINBUF-1) break;	ptr++;

				// length
				LEN1 += 1;
			}

			// end sequence and next
			BUFFER1[LEN1+1] = '\0'; ptr += 1;

			//////////////// check end of buffer /////////////////

			// check read
			if (ptr - BUFFER0 > MAXBUF)
			{
				// correct length
				LEN1 = BUFFER0 - buf + MAXBUF;

				// trunk buffer
				BUFFER1[LEN1] = '\0';

				// update pointer
				ptr = buf + LEN1;
			}

			// check end of buffer
			if (ps2 + ptr - BUFFER0 > sze)
			{
				// correct length
				LEN1 = (sze - ps2) - (buf - BUFFER0);

				// trunk buffer
				BUFFER1[LEN1] = '\0';

				// update pointer
				ptr = buf + LEN1;
			}

			//////////////// read header/sequence ////////////////

			// Check if seq id
			if (BUFFER1[0] == '>')
			{
				// check end of line
				if (strchr(BUFFER1, '\n') == NULL) break;

				// Inc seq
				nb1 += 1;

				// Init len
				ln1 = 0;
			}
			else
			{
				// convert to upper case
				for (__int64 i=0; i<LEN1; i++) SEQ[nb1-1][ln1+i] = upper(BUFFER1[i]); SEQ[nb1-1][ln1+LEN1] = '\0';

				// total length
				nb2 += LEN1;

				// inc length
				ln1 += LEN1;
			}

			// pointer
			buf = ptr;

			// check end of buffer
			if (buf - BUFFER0 >= MAXBUF) break;

			// check end of buffer
			if (ps2 + buf - BUFFER0 >= sze) break;
		}

		// update position
		ps2 += buf - BUFFER0;

		// Verbose
		printf_s("length : %15s\r", format(nb2, ','));
	}

	// Close stream
	fclose(file2);

	// check heap
	if (NBR != nb1) { printf_s("FATAL ERROR : Reading fasta file failed !!!\n"), exit(EXIT_FAILURE); }

	// Verbose
	printf_s("\n");

	// total length
	return nb2;
}

//---------------------------------------------------------------------------
// main routine
void Palindromes(char **SEQ, __int64 *LEN, __int64 NBR, __int64 len, __int64 pos, __int64 seg, __int64 iter)
{
	// Verbose
	printf_s("\n===============================\n");

	// seq counters
	__int64 nbr1 = 0;
	__int64 nbr2 = 0;

	// seg counters
	__int64 nbr3 = 0;
	__int64 nbr4 = 0;

	// pointers
	char *ptr1 = NULL;
	char *ptr2 = NULL;

	// reverse
	char *REV = NULL;

	// private heap
	HANDLE Heap = NULL;

	// suffix tree
	CTree *Tree = NULL;

	//////////////////////////////////////////////

	// Verbose
	printf_s("STEP 3 : Palindromic Sequences\n");

	// Create private heap
	Heap = AllocTree(Tree);

	// Free memory
	for (__int64 i=0; i<NBR; i++)
	{
		// reverse sequence
		REV = _strrev(_strdup(SEQ[i]));

		// complement sequence
		for (__int64 j=0; j<LEN[i]; j++) REV[j] = compl(REV[j]);

		// pointer
		ptr1 = SEQ[i];
		ptr2 = REV + LEN[i] - len;

		// Forward search: add palandromic sequence to suffix tree
		for (__int64 j=0; j<LEN[i]-len+1; j++, nbr1++)
		{
			// check segment first position
			if (++nbr3 < pos) continue;

			// check segment length
			if (++nbr4 > seg) break;

			// Verbose
			if (j % 100000 == 0) printf_s("length : %15s (iter %I64d)\r", format(nbr1, ','), iter);

			// Check filter and add suffix to tree
			SetTree(Tree, ptr1, ptr2, len, 0, j+1, Heap);

			// Next
			ptr1++; ptr2--;
		}

		// verbose
		printf_s("length : %15s (iter %I64d)\r", format(nbr1, ','), iter);

		// free memory
		delete[] REV;

		// check segment length
		if (nbr4 > seg) break;
	}

	//////////////////////////////////////////////
	/////////////////// STEP 4 ///////////////////
	//////////////////////////////////////////////

	// Sort suffix tree
	SortTree(Tree);

	// open output stream
    FILE *file = OpenFile(".\\DATA\\Palindromes.dat", "wb");

	// init buffer
  __int64 MAX = MINBUF;
	char *BUF = new char[MAX];

	// Save suffixes
	GetTree(Tree, 0, nbr2, BUF, MAX, file);

	// Close handle
	fclose(file);

	// free memory
	delete[] BUF;

	// Destroy heap
	HeapDestroy(Heap);

	//////////////////////////////////////////////

	// Verbose
	printf_s("\n===============================\n");
}

//---------------------------------------------------------------------------
// Usage : PALINDROM [file] [start] [length] [palindromes] [executions]
// Arguments :
// [file]        : Fasta format input file
// [start]       : segment start position
// [length]      : segment length
// [palindromes] : palindromes length
// [executions]  : number of executions
void main(int argc, char **argv)
{
	// check args
	if (argc != 6) exit(0);

	//////////////////////////////////////////////
	///////// sequence segment to search /////////
	//////////////////////////////////////////////

	// segment start position
	__int64 pos = _atoi64(argv[2]);

	// segment length in kb
	__int64 seg = _atoi64(argv[3]);

	//////////////////////////////////////////////
	/////// palindromic sequences to search //////
	//////////////////////////////////////////////

	// palindromic length
	__int64 len = _atoi64(argv[4]);

	//////////////////////////////////////////////
	//////////// number of executions ////////////
	//////////////////////////////////////////////

	// number of execution times
	__int64 exe = _atoi64(argv[5]);

	//////////////////////////////////////////////

	// number, total length, sequence lengths, and sequences
	__int64 NBR = 0; __int64 TOT = 0; __int64 *LEN; char **SEQ;

	// Read fasta file format
	ReadFastaFile(argv[1], NBR, LEN, SEQ);

	/////////////////////////////////////////////////

	// init execution time
	__int64 *EXE = new __int64[exe];

	// run N times Palindromic search
	for (__int64 i=0; i<exe; i++)
	{
		// get start time
		auto time1 = std::chrono::high_resolution_clock::now();

		// Search palandromic sequences
		Palindromes(SEQ, LEN, NBR, len, pos, seg*1000, i+1);

		// get last time
		auto time2 = std::chrono::high_resolution_clock::now();

		// get execution time
		EXE[i] = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();
	}

	/////////////////////////////////////////////////

	// init average
	double AVR = 0;

	// init standard deviation
	double STD = 0;

	// calculate average
	for (int i=0; i<exe; i++) AVR += (double)EXE[i]; AVR /= exe;

	// calculate standard deviation
	for (int i=0; i<exe; i++) STD += pow(EXE[i] - AVR, 2.0f); STD = sqrt(STD / exe);

	/////////////////////////////////////////////////

	// Open output file

	FILE *file = OpenFile(".\\DATA\\Times.dat", "wt");

	// separator
	fprintf_s(file, "\n===============================\n");

	// Write input file information
	fprintf_s(file, "Input sequences :\n");
	fprintf_s(file, "Input file          : %s\n", argv[1]);
	fprintf_s(file, "Number of sequences : %I64d\n", NBR);
	fprintf_s(file, "Total length        : %I64d\n", TOT);

	// separator
	fprintf_s(file, "\n===============================\n");

	// Write palindromes information
	fprintf_s(file, "Find palindromes in segment :\n");
	fprintf_s(file, "Palindromes Length  : %I64d\n", len);
	fprintf_s(file, "Segment length      : %I64d\n", seg*1000);
	fprintf_s(file, "Segment start       : %I64d\n", pos);
	fprintf_s(file, "Segment end         : %I64d\n", seg*1000+pos-1);

	// separator
	fprintf_s(file, "\n===============================\n");

	// Write executions information
	fprintf_s(file, "Execution time :\n");
	fprintf_s(file, "Runs : %I64d\n", exe);
	fprintf_s(file, "Mean : %.4lf seconds\n", AVR / 1000000);
	fprintf_s(file, "Stdv : %.4lf seconds\n", STD / 1000000);

	// separator
	fprintf_s(file, "\n===============================\n");

	// cloes stream
	fclose(file);

	/////////////////////////////////////////////////

	// verbose
	printf_s("\n===============================\n");

	// verbose
	printf_s("Check output files \".\\DATA\\Times.dat\" and \".\\DATA\\Palindromes.dat\n\"");

	/////////////////////////////////////////////////

	// Free memory
	for (__int64 i=0; i<NBR; i++)
	delete[] SEQ[i];
	delete[] SEQ;
	delete[] LEN;
	delete[] EXE;
}