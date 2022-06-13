#pragma once
#include "Defines.h"
#include "Control.h"
#include "Geometry.h"
#include "Solver.h"
#include "Xsec.h"

/* Writing output functions */
void PrintOutput();

/* Reading input functions */
void ReadXSLib();
void ReadInput();
void ScanXsecBlock(ifstream& file);
void ScanOptionBlock(ifstream& file);
void ScanGeometryBlock(ifstream& file);
void ScanRadConfCard(ifstream& file, string line);
int GetBlockID(string line);
int GetCardID(int block, string line);

int CountFields(string line, const char* deliminator, bool is_repeat=true);          /// Count the number of divided parts
int SplitFields(string line, string*& fields, const char* deliminator, bool is_repeat=true); /// Split line into several parts
int Repeat(string& field);
void GetLine(ifstream& fin, string& oneline);                                                 /// Util function for IO
void UpperCase(string& line);

int ToInteger(string line);
bool ToBoolean(string line);
double ToDouble(string line);

/* IO variables */
// Blocks and cards
const int num_blocks = 5;
const string BlockName[num_blocks]
{
	"CASEID",
	"XSEC",
	"OPTION",
	"GEOMETRY",
	"."
};

const int num_cards = 8;
const string CardName[num_blocks][num_cards]
{
	/// CASEID
	{},
	/// XSEC
	{ "LIB_TYPE",   "GROUP_SPEC", "DIRECTORY" },
	/// OPTION
	{ "OUTER_CRITERIA", "INNER_CRITERIA", "MESH", "SOLVER", "THREAD"},
	/// GEOMETRY
	{ "PITCH",      "AX_MESH",    "ALBEDO",    "ASSEMBLY",   
	  "RAD_CONF"    },
	/// End Point
	{}
};

/// Enum class
enum class BLOCK{
	CASEID,
	XSEC,
	OPTION,
	GEOMETRY,
};

namespace Card {
	enum class XsecCard {
		LIB_TYPE,
		GROUP_SPEC,
		DIRECTORY,
	};

	enum class OptionCard {
		OUTER_CRITERIA,
		INNER_CRITERIA,
		MESH,
		SOLVER,
		THREAD,
	};

	enum class GeometryCard {
		PITCH,
		AX_MESH,
		ALBEDO,
		ASSEMBLY,
		RAD_CONF,
	};
}