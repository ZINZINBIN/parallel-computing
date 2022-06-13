#include "IO.h"

/* Print output functions */
void PrintOutput()
{
	ofstream file;
	string fn = Control.GetCaseID() + ".out";
	int num_axial = 0, num_vertical = 0, num_horizon = 0;
	int ng = XS.GetNumOfGruop();

	Geometry.GetBasicNum(num_horizon, num_vertical, num_axial);

	file.open(fn, ios::out);
	file.setf(ios_base::floatfield, ios_base::scientific);
	for (int iz = 0; iz < num_axial; iz++) {
		for (int iasyY = 0; iasyY < num_vertical; iasyY++) {
			for (int iasyX = 0; iasyX < num_horizon; iasyX++) {
				int assembly_idx = Geometry.GetAssemblyIdx(iasyX, iasyY);
				if (Geometry.GetAssemblyIdx(iasyX, iasyY) == INVALID) {
					file << INVALID;
					for (int ig = 0; ig < ng; ig++) {
						file << "      " << INVALID;
					}
					file << endl;
				}
				else {
					file << num_horizon * num_vertical * iz + assembly_idx;
					for (int ig = 0; ig < ng; ig++) {
						file << "      " << BiCGSTAB.GetAssemblyFlux(ig, iasyX, iasyY, iz);
					}
					file << endl;
				}
			}
		}
	}
	file.close();
}


/* Read command line to get a input file name */
void Control_t::ReadCMDLine(int argc, char** argv)
{

	// Set input file name
	string fn;
	if (argc < 2) fn = "Input.inp";
	else          fn = argv[1];

	Control.SetInputFile(fn);

	// Get Case ID
	ifstream file;
	string oneline;
	string* fields = NULL;

	int block;

	file.open(fn);
	if (file.fail()) Control.Error(ErrorCode::INPUT_OPEN_ERROR); /// Confirm file existence

	while (!file.eof()) {
		GetLine(file, oneline);
		if (oneline.empty()) continue;
		if (oneline[0] != BLANK) block = GetBlockID(oneline);
		switch (block)
		{
		case (int)BLOCK::CASEID:
			SplitFields(oneline, fields, " ");
			Control.SetCaseID(fields[1]);
			break;

		case (int) BLOCK::OPTION:
			ScanNumThread(file);
			break;
		}
	}

	file.close();
}

// Scan Num Thread
void Control_t::ScanNumThread(ifstream& file)
{
	string oneline;
	string* fields = NULL;
	int block, card;
	int pos;

	while (!file.eof()) {
		pos = (int)file.tellg();
		GetLine(file, oneline);
		Control.Message(oneline, false, false, true);

		if (oneline.empty()) continue;
		block = GetBlockID(oneline);

		if (block != INVALID) {
			file.seekg(pos);
			break;
		}

		card = GetCardID((int)BLOCK::OPTION, oneline);

		switch (card)
		{
		case (int)Card::OptionCard::THREAD:
			SplitFields(oneline, fields, " ");
			_num_omp_threads = ToInteger(fields[1]);
			break;

		}
	}
}

// Read XS Lib
void ReadXSLib()
{
	ifstream file;
	string fn = Control.GetXSFile();
	string oneline;
	string *fields = NULL;

	int ng          = XS.GetNumOfGruop();

	file.open(fn);
	if (file.fail()) Control.Error(ErrorCode::XSLIB_OPEN_ERROR); /// Confirm file existence

	stringstream s;
	s << "XSlib File             : " << Control.GetXSFile();
	Control.Message(s.str());
	Control.Message("Reading XS library...");
	while (!file.eof()) {
		GetLine(file, oneline);
		if (oneline.empty()) continue;
		
		SplitFields(oneline, fields, " ");
		UpperCase(fields[0]);

		if (!fields[0].compare("MATERIAL")) {
			Mat_t new_mat;
			new_mat._index = ToInteger(fields[1]);
			new_mat.SetNumOfGroup(ng);
			/// XS
			for (int ig = 0; ig < ng; ig++) {
				GetLine(file, oneline);
				SplitFields(oneline, fields, " ");

				new_mat._D[ig]      = 1.0 / (3.0 * ToDouble(fields[0]));
				new_mat._sigr[ig]   = ToDouble(fields[1]);
				new_mat._nusigf[ig] = ToDouble(fields[2]);
				new_mat._chi[ig]    = ToDouble(fields[3]);
			}
			/// Scattering XS
			for (int ig = 0; ig < ng; ig++) {
				GetLine(file, oneline);
				SplitFields(oneline, fields, " ");

				for (int igg = 0; igg < ng; igg++) {
					new_mat._scat(ig, igg) = ToDouble(fields[igg]);
				}
			}

			/// Removal
			for (int ig = 0; ig < ng; ig++) {
				for (int igg = 0; igg < ng; igg++) {
					if (igg == ig) continue;
					new_mat._sigr[ig] = new_mat._sigr[ig] + new_mat._scat(ig, igg);
				}
			}

			XS.SetMaterialType(new_mat);
		}
	}

	file.close();
}

// Read input file
void ReadInput()
{
	ifstream file;
	string fn = Control.GetInputFile();
	string oneline;
	string *fields = NULL;

	int block;

	file.open(fn);
	if (file.fail()) Control.Error(ErrorCode::INPUT_OPEN_ERROR); /// Confirm file existence

	Control.Message("Reading input...");
	while (!file.eof()) {
		GetLine(file, oneline);
		Control.Message(oneline, false, false, true);
		if (oneline.empty()) continue;
		if (oneline[0] != BLANK) block = GetBlockID(oneline);
		switch (block)
		{
		case (int) BLOCK::CASEID:
			SplitFields(oneline, fields, " ");
			Control.SetCaseID(fields[1]);
			break;

		case (int) BLOCK::XSEC:
			ScanXsecBlock(file);
			break;

		case (int) BLOCK::OPTION:
			ScanOptionBlock(file);
			break;

		case (int) BLOCK::GEOMETRY:
			ScanGeometryBlock(file);
			break;
		}
	}

	file.close();
}

// Scan XSEC block in an input file
void ScanXsecBlock(ifstream& file)
{
	string oneline;
	string* fields = NULL;
	int block, card;
	int pos;

	while (!file.eof()) {
		pos = (int) file.tellg();
		GetLine(file, oneline);
		Control.Message(oneline, false, false, true);

		if (oneline.empty()) continue;
		block = GetBlockID(oneline);

		if (block != INVALID) {
			file.seekg(pos);
			break;
		}

		card = GetCardID((int)BLOCK::XSEC, oneline);

		switch (card)
		{
		case (int) Card::XsecCard::LIB_TYPE:
			SplitFields(oneline, fields, " ");
			XS.SetLibType(ToInteger(fields[1]));
			break;

		case (int) Card::XsecCard::GROUP_SPEC:
			SplitFields(oneline, fields, " ");
			XS.SetNumGroup(ToInteger(fields[1]));
			break;
		case (int) Card::XsecCard::DIRECTORY:
			SplitFields(oneline, fields, " ");
			Control.SetXSFile(fields[1]);
			break;
		}
	}
}

// Scan OPTION block in an input file
void ScanOptionBlock(ifstream& file)
{
	string oneline;
	string* fields = NULL;
	int block, card;
	int pos;

	while (!file.eof()) {
		pos = (int)file.tellg();
		GetLine(file, oneline);
		Control.Message(oneline, false, false, true);

		if (oneline.empty()) continue;
		block = GetBlockID(oneline);

		if (block != INVALID) {
			file.seekg(pos);
			break;
		}

		card = GetCardID((int)BLOCK::OPTION, oneline);

		switch (card)
		{
		case (int)Card::OptionCard::OUTER_CRITERIA:
			SplitFields(oneline, fields, " ");
			Control.SetEigenCriteria(ToDouble(fields[1]));
			Control.SetPsiCriteria(ToDouble(fields[2]));
			Control.SetMaxOuter(ToInteger(fields[3]));
			break;

		case (int)Card::OptionCard::INNER_CRITERIA:
			SplitFields(oneline, fields, " ");
			Control.SetInnerCriteria(ToDouble(fields[1]));
			Control.SetMaxInner(ToInteger(fields[2]));
			break;

		case (int)Card::OptionCard::MESH: {
			SplitFields(oneline, fields, " ");
			int num_meshX = ToInteger(fields[1]);
			int num_meshY = ToInteger(fields[2]);
			int num_meshZ = ToInteger(fields[3]);
			Geometry.SetNumOfMesh(num_meshX, num_meshY, num_meshZ);
			break;
		}

		case (int) Card::OptionCard::SOLVER:
			SplitFields(oneline, fields, " ");
			UpperCase(fields[1]);
			if (fields[1].compare("BICGSTAB")) Control.SetSolverType(SolverType::BICGSTAB);
			else if (fields[1].compare("SOR")) Control.SetSolverType(SolverType::SOR);
		}
	}
}

// Scan Geometry block in an input file
void ScanGeometryBlock(ifstream& file)
{
	string oneline;
	string *fields = NULL;
	int block, card;
	int pos;
	int num_fields;

	while (!file.eof()) {
		pos = (int) file.tellg();
		GetLine(file, oneline);
		Control.Message(oneline, false, false, true);

		if (oneline.empty()) continue;
		block = GetBlockID(oneline);

		if (block != INVALID) {
			file.seekg(pos);
			break;
		}

		card = GetCardID((int)BLOCK::GEOMETRY, oneline);

		switch (card)
		{
		case (int) Card::GeometryCard::PITCH:
			SplitFields(oneline, fields, " ");
			Geometry.SetAssemblyPitch(ToDouble(fields[1]));
			break;
		
		case (int)Card::GeometryCard::AX_MESH: {
			int num_axial_node = SplitFields(oneline, fields, " ") - 1;
			Geometry.SetNumOfAxialNode(num_axial_node);
			for (int iz = 1; iz < num_axial_node + 1; iz++)
				Geometry.SetAxialNode(iz - 1, ToDouble(fields[iz]));
			break;
		}
		
		case (int) Card::GeometryCard::ALBEDO:
			SplitFields(oneline, fields, " ");
			for (int i = 0; i < 6; i++) Geometry.SetAlbedo(i, ToDouble(fields[1 + i]));
			break;
		
		case (int)Card::GeometryCard::ASSEMBLY: {
			num_fields = SplitFields(oneline, fields, " ");
			Assembly_t new_assembly;
			new_assembly._index = ToInteger(fields[1]);
			new_assembly._axial_material = new int[num_fields - 2];
			for (int iz = 0; iz < num_fields - 2; iz++)
				new_assembly._axial_material[iz] = ToInteger(fields[2 + iz]) - 1;
		
			Geometry.SetAssemblyType(new_assembly);
			break;
		}
		
		case (int) Card::GeometryCard::RAD_CONF:
			ScanRadConfCard(file, oneline);
			break;
		}
	}
}

// Scan XSEC block in an input file
void ScanRadConfCard(ifstream& file, string line)
{
	string oneline;
	string* fields = NULL;
	int block, card;
	int num_fields;
	int pos;

	/// Check radial configuration property
	int angle;
	bool is_cent = false;
	SplitFields(line, fields, " ");
	angle = ToInteger(fields[1]);

	UpperCase(fields[2]);
	if (fields[2].compare("CENT"))      is_cent = true;
	else if (fields[2].compare("EDGE")) is_cent = false;
	Geometry.SetRadConfAngle(angle, is_cent);

	/// Scan radial configuration
	int num_assembly = 0;
	int num_vertical = 0, num_horizon = 0;
	vector<string> onelines;
	while (!file.eof()) {
		pos = (int)file.tellg();
		GetLine(file, oneline);
		Control.Message(oneline, false, false, true);

		if (oneline.empty()) continue;
		block = GetBlockID(oneline);

		if (block != INVALID) {
			file.seekg(pos);
			break;
		}

		card = GetCardID((int)BLOCK::GEOMETRY, oneline);

		if (card != INVALID) {
			file.seekg(pos);
			break;
		}
		num_fields = SplitFields(oneline, fields, " ");

		num_vertical++;
		num_horizon   = max(num_horizon, num_fields);
		num_assembly += num_fields;
		onelines.push_back(oneline);
	}

	int iy = 0, idx = 0;
	Geometry.SetRadConfSize(num_assembly, num_horizon, num_vertical);
	Array_t<int> rad_conf_mat, rad_conf_idx; 
	rad_conf_mat.Create(num_horizon + 2, num_vertical + 2); rad_conf_mat.Offset(-1, -1);
	rad_conf_idx.Create(num_horizon + 2, num_vertical + 2); rad_conf_idx.Offset(-1, -1);
	rad_conf_mat = INVALID; rad_conf_idx = INVALID;
	while (!onelines.empty()) {
		oneline = onelines.front(); onelines.erase(onelines.begin());
		num_fields = SplitFields(oneline, fields, " ");

		int residual = num_horizon - num_fields;

		if (angle == 90) {
			for (int ix = 0; ix < num_fields; ix++) {
				rad_conf_mat(ix, iy) = ToInteger(fields[ix]);
				rad_conf_idx(ix, iy) = idx;
				idx++;
			}
		}
		else if (angle == 360) {
			int half_residual = (int) (((double) residual) / 2.0);
			for (int ix = 0; ix < num_fields; ix++) {
				rad_conf_mat(half_residual + ix, iy) = ToInteger(fields[ix]) - 1;
				rad_conf_idx(half_residual + ix, iy) = idx;
				idx++;
			}
		}
		iy++;
	}
	Geometry.SetRadConfMat(rad_conf_mat);
	Geometry.SetRadConfIdx(rad_conf_idx);
	rad_conf_mat.Destroy();
	rad_conf_idx.Destroy();
}

// Get block name of line
int GetBlockID(string line)
{
	int pos_end  = (int) line.find(BLANK, 0);
	string block = line.substr(0, pos_end);

	UpperCase(block);
	for (int i = 0; i < num_blocks; i++)
		if (!block.compare(BlockName[i])) return i;

	return INVALID;
}

// Get card name of line in a block
int GetCardID(int block, string line)
{
	int pos_start = (int)line.find_first_not_of(BLANK, 0);
	int pos_end = (int) line.find(BLANK, pos_start);
	string card = line.substr(pos_start, pos_end - 1);

	UpperCase(card);
	for (int i = 0; i < num_cards; i++)
		if (!card.compare(CardName[block][i])) return i;

	return INVALID;
}


// Count the number of sub-line divided by deliminator
int CountFields(string line, const char* deliminator, bool is_repeat)
{
	int num_fields = 0;
	int rpt;
	char* context = NULL, *pch;

	pch = strtok((char*)line.c_str(), deliminator, &context); /// Set first part
	if (pch != NULL) {
		while (pch != NULL) {
			string field = pch;
			rpt = 1;
			if (is_repeat) rpt = Repeat(field);
			num_fields += rpt;                                 /// Count the number of fields
			pch = strtok(context, deliminator, &context); /// Find next part
		}
	}
	return num_fields;
}

// Split given string line into sub fields with a given deliminaotr standard
int SplitFields(string line, string*& fields, const char* deliminator, bool is_repeat)
{
	int num_fields = 0;
	char *context = NULL, *pch;

	num_fields = CountFields(line, deliminator, is_repeat); /// Get the number of sub lines  

	if (fields != NULL) delete[] fields;         /// Set fields
	if (!num_fields) return 0;                   /// if incorrect form is given, return 0

	fields = new string[num_fields];

	int i = 0;
	int rpt;
	pch = strtok((char*)line.c_str(), deliminator, &context); /// Find the first part
	while (pch != NULL) {
		string field = pch;
		rpt = 1;
		if (is_repeat) rpt = Repeat(field);
		for (int ir = 0; ir < rpt; ir++) {
			fields[i] = field;
			i++;
		}                                    /// Store a part
		pch = strtok(context, deliminator, &context);         /// Find new part
	}

	return num_fields;
}

// Repeat based on criteria such as '*'
int Repeat(string& field)
{
	int pos = (int) field.find('*');
	if (pos == string::npos) return 1;
	else {
		pos = (int) field.find('*');
		int rpt = ToInteger(field.substr(0, pos));
		field = field.substr(pos + 1, field.size() - pos - 1);
		return rpt;
	}
}

// Reading Factors
void GetLine(ifstream& fin, string& oneline)
{
	getline(fin, oneline);
	replace(oneline.begin(), oneline.end(), '\t', BLANK); ///changing the line to usual form
#ifdef __linux__
	replace(oneline.begin(), oneline.end(), '\r', BLANK); ///changing the line to usual form in linux
#endif
	int pos = (int) oneline.find_first_of(PANG);
	oneline = oneline.substr(0, pos);
	if (oneline.find_first_not_of(BLANK) == string::npos) /// Checking empty
	{
		oneline.clear();
		return; //end of function
	}
}

// Uppercase string
void UpperCase(string& line)
{
	transform(line.begin(), line.end(), line.begin(), ::toupper);
}


// Change string to integer
int ToInteger(string line)
{
	int val;
	try {
		val = stoi(line);
	}
	catch (invalid_argument) {
		Control.Error(ErrorCode::INTEGER_ERROR);
	}

	return val;
}

// Change string to boolean
bool ToBoolean(string line)
{
	UpperCase(line);

	if (line == "T")          return true;
	else if (line == "TRUE")  return true;
	else if (line == "F")     return false;
	else if (line == "FALSE") return false;
	else Control.Error(ErrorCode::BOOLEAN_ERROR);

	return false;
}


// Change string to double 
double ToDouble(string line)
{
	double val;
	try {
		val = stod(line);
	}
	catch (invalid_argument) {
		Control.Error(ErrorCode::DOUBLE_ERROR);
	}

	return val;
}