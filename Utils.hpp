//---------------------------------------------------------------------------
//
#include <chrono>
#include <fstream>
#include <windows.h>

// load fasta file in upper case format
#define lower(c)	((c >= 'A' && c <= 'Z') ? (c - 'A' + 'a') : c)
#define upper(c)	((c >= 'a' && c <= 'z') ? (c - 'a' + 'A') : c)

// sizes
const int MAXBUF = 65536;		// Max output buffer
const int MINBUF = 16384;		// Min buffer length
const int SYMBOL =     4;		// Number of residues

// buffers
char BUFFER0[MAXBUF];			// Output buffer (read/write)
char BUFFER1[MINBUF];			// Output buffer (sequence)
char BUFFER2[MINBUF];			// Output buffer (occurrence)

// buffers test
bool TEST1 = false;				// BUFFER1 test
bool TEST2 = false;				// BUFFER2 test

// buffers length
__int64 LEN1 = 0;				// BUFFER1 length
__int64 LEN2 = 0;				// BUFFER2 length

auto MAXSEQ = 25;				// Max length for long sequences
auto MINSEQ =  4;				// Max length for long sequences

//---------------------------------------------------------------------------
// Swap two variables
template<class typ> inline void SWAP(typ &a, typ &b) { typ t = a; a = b; b = t; }

//---------------------------------------------------------------------------
// Structur : Suffixe Tree
typedef struct CTree
{
	int    NBR;
	char   RES;
  __int64 *POS;
	CTree *BRO;
	CTree *SON;
} CTree;

//---------------------------------------------------------------------------
// Add new to tree
CTree* AddBRO(CTree *Tree, HANDLE heap)
{
	// get first
	CTree *ptr = Tree;

	// get last
	while (ptr->BRO != NULL) { ptr = ptr->BRO; }

	// add new
	ptr->BRO = (CTree*) HeapAlloc(heap, HEAP_ZERO_MEMORY, sizeof(CTree));

	// checck allocation error
	if (ptr->BRO == NULL) { printf_s("MEMORY REALLOCATION EEROR !!!\n"); exit(EXIT_FAILURE); }

	// get new
	ptr = ptr->BRO;
	
	// set no next
	if (ptr != NULL)
	{
		ptr->BRO = NULL;
		ptr->SON = NULL;
	}

	// send new
	return ptr;
}

//---------------------------------------------------------------------------
// Add new to tree
CTree* AddSON(CTree *Tree, HANDLE heap)
{
	// get first
	CTree *ptr = Tree;

	// get last
	while (ptr->SON != NULL) { ptr = ptr->SON; }

	// add new
	ptr->SON = (CTree*) HeapAlloc(heap, HEAP_ZERO_MEMORY, sizeof(CTree));

	// checck allocation error
	if (ptr->SON == NULL) { printf_s("MEMORY REALLOCATION EEROR !!!\n"); exit(EXIT_FAILURE); }

	// get new
	ptr = ptr->SON;
	
	// set no next
	if (ptr != NULL)
	{
		ptr->BRO = NULL;
		ptr->SON = NULL;
	}

	// send new
	return ptr;
}

//---------------------------------------------------------------------------
// memory re-allocation on heap
template<class type1, class type2> inline type1* Realloc(type1 *src, type2 size, HANDLE heap)
{
	// check null
	if (src == NULL)
	{
		// heap memory allocation
		type1 *ptr = (type1*) HeapAlloc(heap, HEAP_ZERO_MEMORY, sizeof(type1) * size);

		// check error
		if (ptr == NULL) { printf_s("MEMORY REALLOCATION EEROR !!!\n"); exit(EXIT_FAILURE); }

		// send new
		return ptr;
	}
	else
	{
		// heap memory reallocation
		type1* ptr = (type1*) HeapReAlloc(heap, HEAP_ZERO_MEMORY, src, sizeof(type1) * size);

		// check error allocation
		if (ptr == NULL) { printf_s("MEMORY REALLOCATION EEROR !!!\n"); exit(EXIT_FAILURE); }

		// send new
		return ptr;
	}
}

//---------------------------------------------------------------------------
// memory re-allocation
template<class type1, class type2> inline type1* Realloc(type1 *src, type2 size1, type2 size2)
{
	// alloc memory
	type1 *dst = new type1[size2];

	// check if null
	if (src == NULL) return dst;

	// copy source
	errno_t err = memcpy_s(dst, size2 * sizeof(type1), src, size1 * sizeof(type1));

	// check error allocation
	if (err != 0) { printf_s("MEMORY REALLOCATION EEROR !!!\n"); exit(EXIT_FAILURE); }

	// free memory
	delete[] src;

	// send new
	return dst;
}

//---------------------------------------------------------------------------
// memory allocation on heap
template<typename type> type* NEW(__int64 size, HANDLE heap)
{
	// Alloc heap header
	type *ptr = (type*) HeapAlloc(heap, HEAP_ZERO_MEMORY, sizeof(type) * size);

	// check error allocation
	if (ptr == NULL) { printf_s("MEMORY REALLOCATION EEROR !!!\n"); exit(EXIT_FAILURE); }

	// new mem
	return ptr;
}

// ---------------------------------------------------------------------------
// Free Suffix Tree
HANDLE AllocTree(CTree *&Tree)
{
	// Create private heap
	HANDLE Heap = HeapCreate(HEAP_GENERATE_EXCEPTIONS | HEAP_NO_SERIALIZE, 0, 0);

	// check heap
	if (Heap == NULL) { printf_s("FATAL ERROR : Memory allocation failed !!!\n"), exit(EXIT_FAILURE); }

	// Alloc heap header
	Tree = (CTree*) HeapAlloc(Heap, HEAP_ZERO_MEMORY, sizeof(CTree));

	// check heap
	if (Tree == NULL) { printf_s("FATAL ERROR : Memory allocation failed !!!\n"), exit(EXIT_FAILURE); }

	// init Suffix Tree
	Tree->NBR   =  0;
	Tree->RES   = '#';
	Tree->BRO   = NULL;
	Tree->SON   = NULL;

	// default
	return Heap;
}

// ---------------------------------------------------------------------------
// convert int64 to char
char *int2char(__int64 val)
{
	// init char
	char str[_MAX_U64TOSTR_BASE2_COUNT];

	// convert int to char
	_i64toa_s(val, str, _countof(str), 10);

	// return char
	return str;
}

//---------------------------------------------------------------------------
// Get str length
__int64 length(char *str)
{
	// pointer
	char *ptr = str;

	while (true)
	{
		switch (*ptr)
		{
			// end of char
			case '\n' : return ptr - str + 1;
			case '\r' : return ptr - str + 1;
			case '\0' : return ptr - str;

			// next char
			default : ptr += 1;
		}
	}
}

//---------------------------------------------------------------------------
// Convert number to text
inline char* format(__int64 value, char sep)
{
	char str1[MINBUF]; static char str2[MINBUF];

	// convert value to text
	if (_i64toa_s(value, str1, MINBUF, 10) != 0) { printf_s("CONVERSION EEROR !!!\n"); exit(EXIT_FAILURE); }

	// scroll text
	for (__int64 i=0, j=0, k=length(str1); i<length(str1); i++, k--)
	{
		// check separation															// update number
		if (k < length(str1) && k % 3 == 0) str2[j++] = sep, str2[j] = '\0';		str2[j++] = str1[i]; str2[j] = '\0';
	}

	return str2;
}

//---------------------------------------------------------------------------
// Get nucleic acid id from its one letter symbol
inline char na2id(char na)
{
	switch (upper(na))
	{
		// standard aa.
		case 'A' : return 0;
		case 'C' : return 1;
		case 'G' : return 2;
		case 'T' : return 3;

		// Invalid residue
		default : return -1;
	}
}

// ---------------------------------------------------------------------------
// Complement na sequence
inline char compl(char na)
{
	switch (upper(na))
	{
		case 'A' : return 'T';
		case 'G' : return 'C';
		case 'C' : return 'G';
		case 'T' : return 'A';
		case 'U' : return 'A';
		case 'R' : return 'Y';
		case 'Y' : return 'R';
		case 'S' : return 'S';
		case 'W' : return 'W';
		case 'K' : return 'M';
		case 'M' : return 'K';
		case 'B' : return 'V';
		case 'V' : return 'B';
		case 'D' : return 'H';
		case 'H' : return 'D';
		case 'N' : return 'N';
		case '-' : return '-';

		// Invalid residue
		default : exit(EXIT_FAILURE);
	}
}

// ---------------------------------------------------------------------------
// Open file stream
FILE *OpenFile(const char *File, const char *p)
{
	// Open for read
	FILE *stream = _fsopen(File, p, _SH_SECURE);
 
	// check if file open
	if(stream == NULL) { printf_s("File \"%s\" not found !!!\n", File); exit(EXIT_FAILURE); }
 
	// Return file
	return stream;
}

// ---------------------------------------------------------------------------
// Add sequence to suffix tree
void SetTree(CTree *Tree, char *SEQ1, char *SEQ2, __int64 LEN, int POS, __int64 COR, HANDLE heap)
{
	// check if palandromic positions
	if (SEQ1[POS] != SEQ2[POS]) return;

	// check if defined nucleotide bases
	if (na2id(SEQ1[POS]) < 0 || na2id(SEQ2[POS]) < 0) return;

	// get first
	CTree *SON = Tree->SON;

	// Scroll over all
	while (SON != NULL)
	{
		// check if match
		if (SON->RES == SEQ1[POS])
		{
			// check length
			if (POS == LEN-1)
			{
				// new position
				SON->POS = Realloc(SON->POS, ++SON->NBR, heap); 

				// add coordinate
				SON->POS[SON->NBR-1] = COR;
			}
			else
			{	// next
				SetTree(SON, SEQ1, SEQ2, LEN, POS + 1, COR, heap);
			}

			// exit
			return;
		}

		// next
		SON = SON->BRO;
	}

	// init new
	CTree *ADD;

	// check if first new
	if (Tree->SON == NULL)
	{
		// add new
		ADD = AddSON(Tree, heap);
	}
	else
	{
		// add next
		ADD = AddBRO(Tree->SON, heap);
	}

	// New
	ADD->NBR = 0;
	ADD->RES = SEQ1[POS];
	ADD->POS = NULL;
	ADD->BRO = NULL;
	ADD->SON = NULL;

	// check length
	if (POS == LEN-1)
	{
		// New
		ADD->NBR = 1;

		// new position
		ADD->POS = NEW<__int64>(1, heap);

		// add coordinate
		ADD->POS[0] = COR;
	}
	else
	{
		SetTree(ADD, SEQ1, SEQ2, LEN, POS + 1, COR, heap);
	}
}

// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree(CTree *Tree, __int64 POS, __int64 &LEN, char *&BUF, __int64 &MAX, FILE *file)
{
	// get first son
	CTree *SON = Tree->SON;

	// Scroll over sons
	while (SON != NULL)
	{
		// Add current residue to motif
		BUFFER1[POS + 0] = SON->RES;
		BUFFER1[POS + 1] = '\0';

		if (SON->NBR > 0)
		{
			// get seq length
			LEN1 = length(BUFFER1);
			LEN2 = 0;

			// init buffer
			sprintf_s(BUF, MAX, "\0");

			// get positions
			for (int i=0; i<SON->NBR; i++)
			{
				// get current position
				char *ptr = int2char(SON->POS[i]);

				// get buffer length
				LEN2 += strlen(ptr);

				// check buffer length
				if (LEN2 + 1 >= MAX) 
				{
					// resize buffer
					BUF = Realloc(BUF, MAX, MAX+1024); MAX += 1024;
				}

				// separator
				if (i > 0)
				{
					// in length
					LEN2 += 1;

					// add separator
					strcat_s(BUF, MAX, ",");
				}

				// position
				strcat_s(BUF, MAX, ptr);
			}

			// checkk output buffer length
			if (LEN1 + LEN2 + LEN + 2 >= MAXBUF)
			{
				// wrtite to file
				fwrite(BUFFER0, sizeof(char), LEN, file);

				// reset buffer
				LEN = 0; BUFFER0[0] = '\0';
			}

			// pointer to buffer
			char *buf = BUFFER0 + LEN;

			// add sequence to buffer
			for (int i=0; i<LEN1; i++) buf[i] = BUFFER1[i]; buf[LEN1] = '\t';

			// add occurrences to buffer
			for (int i=0; i<LEN2; i++) buf[LEN1 + i + 1] = BUF[i]; buf[LEN1 + LEN2 + 1] = '\n'; buf[LEN1 + LEN2 + 2] = '\0';

			// inc buffer
			LEN += LEN1 + LEN2 + 2;
		}

		// Next
		GetTree(SON, POS+1, LEN, BUF, MAX, file);

		// Next
		SON = SON->BRO;
	}

	// Write sequence
	if (Tree->RES == '#' && LEN > 0) fwrite(BUFFER0, sizeof(char), LEN, file);
}

// ---------------------------------------------------------------------------
// Sort suffix tree
void SortTree(CTree *Tree)
{
	// get first
	CTree *SON1 = Tree->SON;

	// Scroll over
	while (SON1 != NULL)
	{
		// get first
		CTree *SON2 = SON1->BRO;

		// Scroll over
		while (SON2 != NULL)
		{
			// check and swap residues
			if (SON1->RES > SON2->RES)
			{
				SWAP(SON1->RES, SON2->RES);
				SWAP(SON1->NBR, SON2->NBR);
				SWAP(SON1->POS, SON2->POS);
				SWAP(SON1->SON, SON2->SON);
			}

			// Next
			SON2 = SON2->BRO;
		}

		// sort next
		SortTree(SON1);

		// Next
		SON1 = SON1->BRO;
	}
}

// ---------------------------------------------------------------------------
// Free Suffix Tree
void EmptyTree(CTree *Tree)
{
	// next
	if (Tree != NULL)
	{
		// Scroll horizontally
		EmptyTree(Tree->BRO);

		// Scroll vertically
		EmptyTree(Tree->SON);

		// free memory
		delete Tree->BRO; delete Tree->SON;
	}
}
