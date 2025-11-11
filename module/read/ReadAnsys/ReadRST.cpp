/* This file is part of Vistle.

   You can use it under the terms of the GNU Lesser General Public License
   version 2.1 or later, see lgpl-2.1.txt.

 * License: LGPL 2+ */

// =============================================================================
// READRFL Klasse zum lesen von ANSYS RFL-Ergebnisfiles (FLOWTRAN)
// -----------------------------------------------------------------------------
// 17.9.2001  BjÃ¶rn Sander
// 16.1.2002  Sergio Leseduarte
// =============================================================================
/* 
const char *dofname[32] = {"UX",      "UY",      "UZ",   "ROTX", "ROTY",    "ROTZ",    "AX",      "AY",
                           "AZ",      "VX",      "VY",   "VZ",   "unused1", "unused2", "unused3", "unused4",
                           "unused5", "unused6", "PRES", "TEMP", "VOLT",    "MAG",     "ENKE",    "ENDS",
                           "EMF",     "CURR",    "SP01", "SP02", "SP03",    "SP04",    "SP05",    "SP06"};

const char *exdofname[28] = {"DENS", "VISC", "EVIS", "COND", "ECON", "LMD1", "LMD2", "LMD3", "LMD4", "LMD5,
                             "LMD6", "EMD1", "EMD2", "EMD3", "EMD4", "EMD5", "EMD6", "PTOT", "TTOT", "PCOE",
                             "MACH", "STRM", "HFLU", "HFLM", "YPLU", "TAUW", "SPHT", "CMUV"};
 */
#include "ReadRST.h"
#include "ANSYS.h"
#include "AnsysConstants.h"
#include <vistle/core/standardattributes.h>
#include <vistle/module/module.h>
#include <vistle/module/reader.h>
#include <vistle/util/byteswap.h>
#include <memory> // added for unique_ptr

// Define byteSwap functions similar to ByteSwap.h
namespace vistle {
template<typename T>
void byteSwap(T &t)
{
    t = vistle::byte_swap<vistle::little_endian, vistle::big_endian, T>(t);
}

template<typename T>
void byteSwap(T *t, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        *t = vistle::byte_swap<vistle::little_endian, vistle::big_endian, T>(*t);
        ++t;
    }
}
} // namespace vistle
extern const char *dofname[];
extern const char *exdofname[];
//#include <util/coviseCompat.h>

// =============================================================================
// Konstruktor / Destruktor
// =============================================================================

// -----------------------------------------------
// Konstruktor
// -----------------------------------------------

//using namespace vistle;
const double ReadRST::DImpossible_ = pow(2.0, 100.0);
const float ReadRST::FImpossible_ = (float)pow(2.0, 100.0);

/* ReadRST::ReadRST(void) //TODO: evaluate if default constructor can be used instead
{
    DOFData_ = NULL;
    DerivedData_ = NULL;
    rfp_ = NULL;
    ety_ = NULL;
    node_ = NULL;
    element_ = NULL;
    anznodes_ = 0;
    anzety_ = 0;
    anzelem_ = 0;
    nodeindex_ = NULL;
    elemindex_ = NULL;
    timetable_ = NULL;
    SwitchEndian_ = DO_NOT_SWITCH;
    mmap_flag_ = 0;
    mode64_ = false;
} */

// -----------------------------------------------
// Destruktor
// -----------------------------------------------
/* ReadRST::~ReadRST()
{
    Reset(RADIKAL);
} */

// =============================================================================
// Interne Methoden
// =============================================================================

// -----------------------------------------------
// Reset
// -----------------------------------------------
int ReadRST::Reset(int message)
{
    if (!(message & PRESERVE_DOF_DATA)) {
        delete DOFData_;
        DOFData_ = NULL;
    }

    delete DerivedData_;
    DerivedData_ = NULL;

    if (message & RADIKAL) {
        delete DOFData_;
        DOFData_ = NULL;

        delete[] node_;
        node_ = NULL;
        anznodes_ = 0;

        ety_.clear();
        anzety_ = 0;

        element_.clear();
        anzelem_ = 0;

        memset(&header_, 0, sizeof(BinHeader));
        memset(&rstheader_, 0, sizeof(RstHeader));

        if (rfp_ != NULL) {
            fclose(rfp_);
            rfp_ = NULL;
        }
        delete[] nodeindex_;
        nodeindex_ = NULL;
        delete[] elemindex_;
        elemindex_ = NULL;
        delete[] timetable_;
        timetable_ = NULL;

        solheader_.clean();
    }

    return (0);
}

// -----------------------------------------------
// OpenFile
// -----------------------------------------------
// Oeffnet den Ergebnisfile, liest den header und setzt
// den Read-File-Pointer rfp;
// Rueckgabe:
// 0    : alles OK
// 1    : File not found
// 2    : could not read header
// 3    : Read Error Nodal equivalence
// 4    : Read Error Element equivalence
// 5    : Read Error Time table
int ReadRST::OpenFile(const std::string &filename)
{
    if (rfp_ && nodeindex_ && rstheader_.numsets_ != 0)
        return 0;

    if ((rfp_ = fopen(filename.c_str(), "rb")) == NULL) {
        return (1);
    }

    file_des_ = fileno(rfp_);

    if (file_des_ == -1) {
        std::cerr << "fileno failed" << std::endl;
        return (1);
    }

    if (get_file_size() != 0) {
        std::cerr << "Could not get file size" << std::endl;
        return (1);
    }

    std::cout << "file size: " << file_size_ << std::endl;

    // File offen und OK
    // erst mal zwei INTs lesen (lead-in):
    // int 1 ergibt die Laenge des folgenden Records in bytes
    // int 2 ergibt den Typ des Records (z.B. 0x64==Integer)

    /****************************************************************************/
    /* Notes:                                                                                                           */
    /*                                                                                                                      */
    /* ANSYS binary files are written in Fortran which means that data is written   */
    /* blockwise.  Each block or record is delimited by markers, namely a header  */
    /* and trailer, both of them being identical and usually 4 bytes long.  These     */
    /* markers indicate the size of the record.  ANSYS uses particularly variable-  */
    /* length records, however it has been observed that the format of the          */
    /* binary files varies depending on the ANSYS version.                                     */
    /*                                                                                                                      */
    /* In older versions of ANSYS (i.e. V5.6) the size of the records is denoted      */
    /* in BYTES while the second integer after the record header corresponds to  */
    /* the size of the Standard ANSYS File Header which is 100 or 0x64 integers  */
    /* long.  This value is followed by the actual 100 values (integers) or items      */
    /* contained in the ANSYS File Header.                                                             */
    /*                                                                                                                     */
    /* In newer versions (i.e. V11.0) it has been observed that the size of the      */
    /* records is expressed in INTs, that is 0x64 items.  The second integer has  */
    /* the value 2147483648 or 0x80000000 which may be a reference to the    */
    /* data type of the elements contained in the ANSYS File Header, namely        */
    /* integers (signed int ?).  Nevertheless it seems this integer is not counted   */
    /* in the record size.  For more information refer to "Programmer's                */
    /* Manual for Mechanical APDL" (ANSYS Release 12.1)                                     */
    /***************************************************************************/

    // allocate header buffer with RAII
    auto header_buf = std::make_unique<int32_t[]>(100);
    char version[150];

    // Reading the first 4 bytes of the binary file, i.e. the record header: size of the record
    // It is followed by -2147483648, the most negative number that can be expressed by unsigned ints
    /* int32_t num_items;
    fread(&num_items, sizeof(int32_t), 1,
          rfp_); */ // ATTENTION: if you erase this line, you have to add an offset of 1 to the buffer buf below

    int expected_size = 100;
    /* if (num_items != expected_size) {
        std::cerr << "Error reading header. Is not of expected size: " << expected_size << std::endl;
        return (2);
    }

    int header_offset = 0;

    // Reading the header
    if (fread(header_buf.get() + header_offset, sizeof(int32_t), 100, rfp_) != 100) {
        return 2;
    } */
    ReadFortranRecordInts(header_buf.get(), 100, 1);

    std::cout << "First record header value in 100 buf: " << header_buf[0] << std::endl;

    int subversion;
    version[4] = 0;
    memcpy(version, &header_buf[10], 4);
    SwitchEndian(version, 4);

    int ret = sscanf(version, "%d.%d", &header_.version_, &subversion);

    if (ret != 2) {
        header_.version_ = 9;
    }
    std::cout << "File generated by ANSYS V" << header_.version_ << "." << subversion << std::endl;

    if (header_.version_ != 25) {
        // in version 11 and later the record size is given in INTs, not in BYTES
        std::cout << "WARNING: Version may not be supported. Code only tested for ANSYS V25." << std::endl;
    }

    // assign header fields...
    header_.filenum_ = SwitchEndian(header_buf[1]);
    header_.format_ = SwitchEndian(header_buf[2]);
    header_.time_ = SwitchEndian(header_buf[3]);
    header_.date_ = SwitchEndian(header_buf[4]);
    header_.unit_ = SwitchEndian(header_buf[5]);
    header_.ansysdate_ = SwitchEndian(header_buf[11]);
    memcpy(header_.machine_, &header_buf[12], 3 * sizeof(int));
    header_.machine_[12] = '\0';
    SwitchEndian(header_.machine_, 12);
    memcpy(header_.jobname_, &header_buf[31], 8 * sizeof(int));
    header_.jobname_[8] = '\0';
    SwitchEndian(header_.jobname_, 8);
    memcpy(header_.product_, &header_buf[17], 2 * sizeof(int));
    header_.product_[8] = '\0';
    SwitchEndian(header_.product_, 8);
    memcpy(header_.label_, &header_buf[19], 1 * sizeof(int));
    header_.label_[4] = '\0';
    SwitchEndian(header_.label_, 4);
    memcpy(header_.user_, &header_buf[20], 3 * sizeof(int));
    header_.user_[12] = '\0';
    SwitchEndian(header_.user_, 12);
    memcpy(header_.machine2_, &header_buf[23], 3 * sizeof(int));
    header_.machine2_[12] = '\0';
    SwitchEndian(header_.machine2_, 12);
    header_.recordsize_ = SwitchEndian(header_buf[26]);
    header_.maxfilelen_ = SwitchEndian(header_buf[27]);
    header_.maxrecnum_ = SwitchEndian(header_buf[28]);
    memcpy(header_.title_, &header_buf[41], 20 * sizeof(int));
    header_.title_[80] = '\0';
    SwitchEndian(header_.title_, 80);
    memcpy(header_.subtitle_, &header_buf[61], 20 * sizeof(int));
    header_.subtitle_[80] = '\0';
    SwitchEndian(header_.subtitle_, 80);

    std::cout << "ANSYS File Header filenumber: " << header_.filenum_ << std::endl;
    std::cout << "File format: " << header_.format_ << std::endl;
    std::cout << "time in compact form: " << header_.time_ << std::endl;
    std::cout << "date in compact form: " << header_.date_ << std::endl;
    std::cout << "Units of measurement: " << header_.unit_ << std::endl;
    std::cout << "ANSYS date in compact form: " << header_.ansysdate_ << std::endl;
    std::cout << "User: " << header_.user_ << std::endl;
    std::cout << "Machine2: " << header_.machine2_ << std::endl;
    std::cout << "Record size: " << header_.recordsize_ << std::endl;
    std::cout << "Max file length: " << header_.maxfilelen_ << std::endl;
    std::cout << "Max record number: " << header_.maxrecnum_ << std::endl;
    std::cout << "Machine: " << header_.machine_ << "   Jobname: " << header_.jobname_ << std::endl;
    std::cout << "Product: " << header_.product_ << "   Label: " << header_.label_ << std::endl;
    std::cout << "Title: " << header_.title_ << "   Subtitle: " << header_.subtitle_ << std::endl;

    // read RST-File header
    int head_size = 3; // number of lead/trailer ints your file uses
    size_t rsthead_size = 80; // number of ints in RST header
    auto buf_rst = std::make_unique<int32_t[]>(rsthead_size);
    if (ReadFortranRecordInts(buf_rst.get(), rsthead_size, head_size) != rsthead_size) {
        std::cout << "error reading rfp_ rst file header" << std::endl;
        return 2;
    }
    //auto buf_rst = std::make_unique<int32_t[]>(73);
    //memcpy(buf_rst.get(), buf_rst_vec.data(), 73 * sizeof(int32_t));

    std::cout << "raw buf[0]=" << buf_rst[0] << " switched=" << SwitchEndian(buf_rst[0]) << std::endl;

    // fun12 - unit number (result file is 12)
    int isResultFile = 12;
    if (header_.filenum_ != isResultFile) {
        ChangeSwitch();
        if (header_.filenum_ == isResultFile) {
            // file can be read using byteswap, prepare to try again
            fclose(rfp_);
            std::cout << "File is byteswapped, trying again ..." << std::endl;
            return OpenFile(filename);
        } else {
            std::cerr << "Cannot correctly read file header. Seems not to be a result file." << std::endl;
            return 2;
        }
    }

    // assign rstheader_ fields ...
    rstheader_.fun12_ = SwitchEndian(buf_rst[1]);
    rstheader_.maxnodes_ = SwitchEndian(buf_rst[2]);
    rstheader_.usednodes_ = SwitchEndian(buf_rst[3]);
    rstheader_.maxres_ = SwitchEndian(buf_rst[4]);
    rstheader_.numdofs_ = SwitchEndian(buf_rst[5]);
    rstheader_.maxelement_ = SwitchEndian(buf_rst[6]);
    rstheader_.numelement_ = SwitchEndian(buf_rst[7]);
    rstheader_.analysis_ = SwitchEndian(buf_rst[8]);
    rstheader_.numsets_ = SwitchEndian(buf_rst[9]);
    rstheader_.ptr_eof_ = SwitchEndian(buf_rst[10]);
    rstheader_.ptr_dsi_ = SwitchEndian(buf_rst[11]);
    rstheader_.ptr_time_ = SwitchEndian(buf_rst[12]);
    rstheader_.ptr_load_ = SwitchEndian(buf_rst[13]);
    rstheader_.ptr_elm_ = SwitchEndian(buf_rst[14]);
    rstheader_.ptr_node_ = SwitchEndian(buf_rst[15]);
    rstheader_.ptr_geo_ = SwitchEndian(buf_rst[16]);
    rstheader_.units_ = SwitchEndian(buf_rst[20]);
    rstheader_.numsectors_ = SwitchEndian(buf_rst[21]);
    rstheader_.ptr_end_ = 0;

    if (SwitchEndian_ == DO_NOT_SWITCH) {
        rstheader_.ptr_end_ = (long long)(buf_rst[23]) << 32 | buf_rst[24];
    } else {
        rstheader_.ptr_end_ = (long long)(SwitchEndian(buf_rst[24])) << 32 | SwitchEndian(buf_rst[23]);
    }

    std::cout << "RST File Header fun12: " << rstheader_.fun12_ << std::endl;
    std::cout << "Max nodes: " << rstheader_.maxnodes_ << std::endl;
    std::cout << "Used nodes: " << rstheader_.usednodes_ << std::endl;
    std::cout << "Max results: " << rstheader_.maxres_ << std::endl;
    std::cout << "Number of DOFs: " << rstheader_.numdofs_ << std::endl;
    std::cout << "Max elements: " << rstheader_.maxelement_ << std::endl;
    std::cout << "Number of elements: " << rstheader_.numelement_ << std::endl;
    std::cout << "Analysis type: " << rstheader_.analysis_ << std::endl;
    std::cout << "Number of result sets: " << rstheader_.numsets_ << std::endl;
    std::cout << "Units of measurement: " << rstheader_.units_ << std::endl;
    std::cout << "Number of sectors: " << rstheader_.numsectors_ << std::endl;

    //understand pointer:
    std::cout << "Pointer to end of file: " << rstheader_.ptr_end_ << std::endl;
    std::cout << "Pointer to Data Step Index Table: " << rstheader_.ptr_dsi_ << std::endl;
    std::cout << "Pointer to Time Table: " << rstheader_.ptr_time_ << std::endl;
    std::cout << "Pointer to Load Step Table: " << rstheader_.ptr_load_ << std::endl;
    std::cout << "Pointer to Element Equivalence Table: " << rstheader_.ptr_elm_ << std::endl;
    std::cout << "Pointer to Nodal Equivalence Table: " << rstheader_.ptr_node_ << std::endl;
    std::cout << "Pointer to Geometry Description: " << rstheader_.ptr_geo_ << std::endl;
    // Header fertig gelesen

    // Jetzt noch die Indextabellen laden
    int true_used_nodes = rstheader_.usednodes_;

    int offset = (rstheader_.ptr_node_ + PTR_OFFSET) * sizeof(int);
    fseek(rfp_, offset, SEEK_SET); // go to node index table position in stream

    nodeindex_ = new int[true_used_nodes];
    int expectedInt = ReadFortranRecordInts(nodeindex_, true_used_nodes, 0);
    if (expectedInt != true_used_nodes) {
        return (3);
    }
    std::cout << "First 10 entries of Node Equivalence Table:" << std::endl;
    for (int i = 0; i < 20 && i < true_used_nodes; ++i) {
        std::cout << "Node " << i << ": " << nodeindex_[i] << std::endl;
    }

    // Elementindex-Tabelle lesen
    elemindex_ = new int[rstheader_.numelement_];
    if (ReadFortranRecordInts(elemindex_, rstheader_.numelement_, head_size) != rstheader_.numelement_) {
        return (3);
    }
    // print first elements from element index table:
    std::cout << "First 10 entries of Element Equivalence Table:" << std::endl;
    for (int i = 0; i < 10 && i < rstheader_.numelement_; ++i) {
        std::cout << "Element " << i << ": " << elemindex_[i] << std::endl;
    }

    // Timetable lesen
    timetable_ = new double[rstheader_.maxres_];
    int time_offset = (rstheader_.ptr_time_ + PTR_OFFSET) * sizeof(int);
    fseek(rfp_, time_offset, SEEK_SET); // go to time table position in stream
    if (DoubleRecord(timetable_, rstheader_.maxres_) != rstheader_.maxres_) {
        return (6);
    }
    std::cout << "First 10 entries of Time Table:" << std::endl;
    for (int i = 0; i < 10 && i < rstheader_.maxres_; ++i) {
        std::cout << "Time " << i << ": " << timetable_[i] << std::endl;
    }

    // That's all folks
    return (0);
}

int ReadRST::IntRecord(int *buf, int len)
{
    int ret = len;
    if (mmap_flag_) {
        memcpy(buf, (char *)mmap_ini_ + actual_off_, len * sizeof(int));
    } else {
        ret = (int)fread(buf, sizeof(int), len, rfp_);
        if (ret != len)
            return ret;
    }
    int item;
    for (item = 0; item < len; ++item) {
        buf[item] = SwitchEndian(buf[item]);
    }
    return ret;
}

// TODO: leads to segmentation fault
int ReadRST::DoubleRecord(double *buf, int len)
{
    int ret = len;
    if (mmap_flag_) {
        memcpy(buf, (char *)mmap_ini_ + actual_off_, len * sizeof(double));
    } else {
        ret = (int)fread(buf, sizeof(double), len, rfp_);
        if (ret != len)
            return ret;
    }
    int item;
    for (item = 0; item < len; ++item) {
        buf[item] = SwitchEndian(buf[item]);
    }
    return ret;
}

// -----------------------------------------------
// SwitchEndian
// -----------------------------------------------
int ReadRST::SwitchEndian(int value)
{
    if (SwitchEndian_ == SWITCH) {
        int ret = 0;
        int bytes[4];

        bytes[0] = (value & 0x000000FF) << 24;
        bytes[1] = (value & 0x0000FF00) << 8;
        bytes[2] = (value & 0x00FF0000) >> 8;
        bytes[3] = (value & 0xFF000000) >> 24;

        ret = bytes[0] | bytes[1] | bytes[2] | bytes[3];
        return (ret);
    }
    return (value);
}

unsigned int ReadRST::SwitchEndian(unsigned int value)
{
    if (SwitchEndian_ == SWITCH) {
        int ret = 0;
        int bytes[4];

        bytes[0] = (value & 0x000000FF) << 24;
        bytes[1] = (value & 0x0000FF00) << 8;
        bytes[2] = (value & 0x00FF0000) >> 8;
        bytes[3] = (value & 0xFF000000) >> 24;

        ret = bytes[0] | bytes[1] | bytes[2] | bytes[3];
        return (ret);
    }
    return (value);
}

void ReadRST::SwitchEndian(char *buf, int length)
{
    if (SwitchEndian_ == DO_NOT_SWITCH) {
        int index(0);
        while (index < length) {
            char tmp1 = buf[index];
            char tmp2 = buf[index + 1];
            buf[index] = buf[index + 3];
            buf[index + 1] = buf[index + 2];
            buf[index + 2] = tmp2;
            buf[index + 3] = tmp1;
            index += 4;
        }
    }
}

void ReadRST::ChangeSwitch()
{
    if (SwitchEndian_ == SWITCH) {
        SwitchEndian_ = DO_NOT_SWITCH;
    } else {
        SwitchEndian_ = SWITCH;
    }
}

double ReadRST::SwitchEndian(double dval)
{
    if (SwitchEndian_ == SWITCH) {
        vistle::byteSwap(dval);
    }
    return dval;
}

// =============================================================================
// externe Methoden
// =============================================================================

// -----------------------------------------------
// Read
// -----------------------------------------------
int ReadRST::Read(const std::string &filename, int num, std::vector<int> &codes)
{
    Reset(0); // alles lÃ¶schen, wenn nÃ¶tig
    // Mal keine Fehlerabfrage, FIXME!
    int problems = OpenFile(filename);
    switch (problems) {
    // 0    : alles OK
    // 1    : File not found
    // 2    : could not read header
    // 3    : Read Error Nodal equivalence
    // 4    : Read Error Element equivalence
    // 5    : Read Error Time table
    case 0:
        //      printf("Open file OK\n");
        break;

    case 1:
        //sendError("Open file: file not found \n");
        std::cerr << "Open file: file not found" << std::endl;
        break;

    case 2:
        std::cerr << "Open file: Read Error, header" << std::endl;
        break;

    case 3:
        std::cerr << "Open file: Read Error, nodal equ tabular" << std::endl;
        break;

    case 4:
        std::cerr << "Open file: Read Error, element equi tabular" << std::endl;
        break;

    case 5:
        std::cerr << "Open file: Read Error, time table" << std::endl;
        break;
    }
    if (problems)
        return problems;
    problems = GetDataset(num, codes);
    switch (problems) {
    // 1        : File ist nicht offen/initialisiert
    // 2        : Read Error DSI-Tabelle
    // 3        : Num ist nicht im Datensatz
    // 4        : Read Error Solution Header
    // 5        : Read Error DOFs
    // 6        : Read Error exDOFs
    case 0:
        //      printf("GetData : OK\n");
        break;

    case 1:
        std::cerr << "GetData : file not open" << std::endl;
        break;

    case 2:
        std::cerr << "GetData : read error: DSI" << std::endl;
        break;

    case 3:
        std::cerr << "GetData : num with value " << num << " exceeds limits" << std::endl;
        break;

    case 4:
        std::cerr << "GetData : read error solution header" << std::endl;
        break;

    case 5:
        std::cerr << "GetData : read error DOFs" << std::endl;
        break;

    case 6:
        std::cerr << "GetData : read error exDOF" << std::endl;
        break;
    }
    if (problems)
        return problems;
    problems = GetNodes();
    std::cout << "GetNodes finished with problems code: " << problems << std::endl;
    switch (problems) {
    // 1        : Read Error Geometrieheader
    // 2        : Read Error Nodes
    // 3        : Read Error Elementbeschreibung
    // 4        : Read Error ETYs
    // 5        : Read Error Elementtabelle
    // 6        : Read Error Elemente
    case 0:
        //      printf("GetNodes : ok\n");
        break;

    case 1:
        std::cerr << "GetNodes : read error geo" << std::endl;
        break;

    case 2:
        std::cerr << "GetNodes : read error nodes" << std::endl;
        break;

    case 3:
        std::cerr << "GetNodes : read error element description" << std::endl;
        break;

    case 4:
        std::cerr << "GetNodes : read error ety" << std::endl;
        break;

    case 5:
        std::cerr << "GetNodes : read error element tabular" << std::endl;
        break;

    case 6:
        std::cerr << "GetNodes : read error elements" << std::endl;
        break;
    }
    //  fclose(rfp);
    //  rfp=NULL;
    return (problems);
}

// -----------------------------------------------
// GetDataset
// -----------------------------------------------
// Liest die Ergebnisdaten fuer einen Datensatz aus dem File
// Rueckgabe:
// 0        : alles OK
// 1        : File ist nicht offen/initialisiert
// 2        : Read Error DSI-Tabelle
// 3        : Num ist nicht im Datensatz
// 4        : Read Error Solution Header
// 5        : Read Error DOFs
// 6        : Read Error exDOFs
int ReadRST::GetDataset(int num, std::vector<int> &codes)
{
    //  int *buf=NULL;
    int size, i, j, sumdof;
    long long offset;
    std::unique_ptr<double[]> dof;

    // File sollte offen sein
    if (rfp_ == NULL)
        return (1);

    // Out of range check
    if (num > rstheader_.numsets_ /*rstheader_.maxres_*/) {
        std::cout << "ReadSHDR: num " << num << " exceeds limit of " << rstheader_.numsets_ << std::endl;
        return (3);
    }

    // Eventuell alte Daten LÃ¶schen
    // DOF-Liste loeschen
    delete DOFData_;
    DOFData_ = NULL;
    /* FIXME ??????????????????? !!!!!!!!!!!!!!
     delete [] node_;
     node_ = NULL;
     anznodes_=0;
     delete [] ety_;
     ety_ = NULL;
     anzety_=0;
     delete [] element_;
     element_ = NULL;
     anzelem_=0;
   */
    // so, alles gelÃ¶scht
    ReadSHDR(num);

    // jetzt Daten einlesen
    // but first read record length...

    if (solheader_.ptr_nodalsol_ == 0) // no DOF data at all
    {
        solheader_.numnodesdata_ = 0;
        DOFData_ = new DOFData; //CreateNewDOFList();
        DOFData_->anz_ = 0;
        DOFData_->nodesdataanz_ = 0;
        //    DOFData_->data_ = new double[0];
        //    DOFData_->nodesdata_ = new int[0];
        return 0;
    }

    offset = solheader_.offset_ + solheader_.ptr_nodalsol_ * 4;
    fseek(rfp_, (long)offset, SEEK_SET);
    int front[2];
    if (IntRecord(front, 2) != 2) {
        return (5);
    }
    sumdof = solheader_.numdofs_ + solheader_.numexdofs_;
    if (header_.version_ < 10) {
        solheader_.numnodesdata_ = (front[0] - 4) / (sizeof(double) * sumdof);
    } else {
        solheader_.numnodesdata_ = solheader_.numnodes_;
    }

    /*
   cerr << "SHD numnodes "<<solheader_.numnodes_ << ' '
        << "SHD numnodesdata "<<solheader_.numnodesdata_ << endl;
   */

    int frontTeoric = sizeof(int) + sizeof(double) * solheader_.numnodes_ * sumdof;

    size = solheader_.numnodesdata_ * (sumdof);
    dof = std::make_unique<double[]>(size);
    // TODO: this DoubleRecord leads to segmentation fault, find out why
    if (DoubleRecord(dof.get(), size) != size) {
        return (5);
    }

    // Erst mal pro DOF einen double-Array erstellen
    DOFData_ = new DOFData; //CreateNewDOFList();
    DOFData_->anz_ = (int)codes.size() * solheader_.numnodesdata_;
    DOFData_->data_ = new double[DOFData_->anz_];
    DOFData_->displacements_ = new double[3 * solheader_.numnodesdata_];
    DOFData_->nodesdataanz_ = solheader_.numnodesdata_;
    DOFData_->nodesdata_ = new int[DOFData_->nodesdataanz_];

    // Be careful, may be there is no output for all nodes!!!
    if (header_.version_ > 9) {
        front[0] = frontTeoric;
    } else if (front[0] < frontTeoric) // FIXME !!!!!!!!!!!!!!!!!!!
    {
        fseek(rfp_, 3 * sizeof(int), SEEK_CUR); // jump over tail of the last record
        if (IntRecord(DOFData_->nodesdata_, DOFData_->nodesdataanz_) != DOFData_->nodesdataanz_) {
            return 5;
        }
    }

    // fill displacements
    for (i = 0; i < 3; ++i) {
        int code_order;
        for (code_order = 0; code_order < solheader_.numdofs_; ++code_order) {
            if (solheader_.dof_[code_order] == (i + 1))
                break;
        }
        if (!(solheader_.mask_ &
              0x200)) // Nur Teilbereich der Knoten wurde verwendet, assume all nodes have output -> FIXME
        {
            int base = i * solheader_.numnodesdata_;
            if (code_order == solheader_.numdofs_) {
                for (j = 0; j < solheader_.numnodesdata_; ++j) {
                    DOFData_->displacements_[base + j] = 0.0;
                }
            } else {
                for (j = 0; j < solheader_.numnodesdata_; ++j) {
                    DOFData_->displacements_[base + j] = dof[j * sumdof + code_order];
                }
            }
        } else {
            int datanum;
            int base = i * solheader_.numnodesdata_;
            for (datanum = 0; datanum < solheader_.numnodesdata_; ++datanum) {
                DOFData_->displacements_[datanum + base] = 0.0;
            }
        }
    }

    // fill requested DOFs
    for (i = 0; i < codes.size(); ++i) {
        DOFData_->dataset_ = num;
        int code_order;
        if (codes[i] < vistle::ansys::EX_OFFSET) // non-extra normal scalar DOF
        {
            // find position for this code
            for (code_order = 0; code_order < solheader_.numdofs_; ++code_order) {
                if (solheader_.dof_[code_order] == codes[i])
                    break;
            }
            if (code_order == solheader_.numdofs_) {
                std::cerr << "The code of a degree of freedom was not found" << std::endl;
                return 5;
            }

            DOFData_->typ_ = codes[i];
            DOFData_->exdof_ = false;
        } else {
            // extra DOF
            DOFData_->typ_ = codes[i] - vistle::ansys::EX_OFFSET;
            for (code_order = 0; code_order < solheader_.numexdofs_; ++code_order) {
                if (solheader_.exdof_[code_order] == DOFData_->typ_)
                    break;
            }
            // !!!!!!!
            if (code_order == solheader_.numdofs_)
                abort();
            code_order += solheader_.numdofs_; // ??????? !!!!!!!!
            DOFData_->exdof_ = true;
        }
        //    cout << "Nodes !solheader_.mask_: "<< !(solheader_.mask_ & 0x200) << endl;
        if (!(solheader_.mask_ &
              0x200)) // Nur Teilbereich der Knoten wurde verwendet, assume all nodes have output -> FIXME
        {
            int base = i * solheader_.numnodesdata_;
            for (j = 0; j < solheader_.numnodesdata_; ++j) {
                DOFData_->data_[base + j] = dof[j * sumdof + code_order];
            }
        } else {
            int datanum;
            for (datanum = 0; datanum < DOFData_->anz_; ++datanum) {
                DOFData_->data_[datanum] = DImpossible_;
            }
        }
        // ACHTUNG: Extradofs werden falsch bezeichnet!
        // ACHTUNG: Bei Teilbereich sind alle Daten Null!
    }
    // dof freed automatically
    // Alle dofs gelesen und gespeichert!

    return (0);
}

// -----------------------------------------------
// ReadSHDR
// -----------------------------------------------
int ReadRST::ReadSHDR(int num)
{
    unsigned int solbuf[103];
    int size, i;
    long long offset;
    double dsol[100];
    std::unique_ptr<unsigned int[]> buf;

    //  SOLUTIONHEADER shdr;

    // File sollte offen sein
    if (rfp_ == NULL)
        return (1);

    // Out of range check
    if (num > rstheader_.numsets_ /* rstheader_.maxres_ */)
        std::cout << "ReadSHDR: num " << num << " exceeds limits " << rstheader_.numsets_ << std::endl;
    return (3);

    // Springe erst mal zu der DSI-Tabelle
    offset = rstheader_.ptr_dsi_ * 4;
    fseek(rfp_, (long)offset, SEEK_SET); // im pointer steht die Anzahl der int-Elemente vom Anfang

    // Jetzt sollte man die Tabelle Lesen koennen. Die ist mitunter aber recht gross
    // gewoehnlich beinhaltet sie 2*1000 Eintraege (besser 1000 64-bit Pointer)
    // dazu kommt dann immer noch der Lead-in (2 Ints) und er Lead-out (1 int)
    size = 2 * rstheader_.maxres_ + 3;
    std::unique_ptr<int[]> buf_tab = std::make_unique<int[]>(size);

    if (fread(buf_tab.get(), sizeof(int), size, rfp_) != size) {
        return (2);
    }
    // jetzt mal diese Tabelle auswerten
    solheader_.offset_ = 0;
    // Hi/Lo lesen umdrehen und einfuegen
    if (SwitchEndian_ == DO_NOT_SWITCH) {
        solheader_.offset_ = ((long long)(buf_tab[num + 2 + rstheader_.maxres_]) << 32 | buf_tab[num + 1]) * 4;
        if (num < rstheader_.numsets_) {
            solheader_.next_offset_ = ((long long)(buf_tab[num + 3 + rstheader_.maxres_]) << 32 | buf_tab[num + 2]) * 4;
        } else {
            solheader_.next_offset_ = file_size_;
        }
    } else {
        solheader_.offset_ =
            ((long long)(SwitchEndian(buf_tab[num + 2 + rstheader_.maxres_])) << 32 | SwitchEndian(buf_tab[num + 1])) *
            4;
        if (num < rstheader_.numsets_) {
            solheader_.next_offset_ = ((long long)(SwitchEndian(buf_tab[num + 3 + rstheader_.maxres_])) << 32 |
                                       SwitchEndian(buf_tab[num + 2])) *
                                      4;
        } else {
            solheader_.next_offset_ = file_size_;
        }
    }
    // jetzt da hin springen und dort einen Solution Header lesen
    fseek(rfp_, (long)solheader_.offset_, SEEK_SET);
    if (fread(solbuf, sizeof(int), 103, rfp_) != 103) {
        return (4);
    }
    // Jetzt die Werte decodieren und zuweisen
    solheader_.numelements_ = SwitchEndian(solbuf[3]);
    solheader_.numnodes_ = SwitchEndian(solbuf[4]);
    solheader_.mask_ = SwitchEndian(solbuf[5]);
    solheader_.loadstep_ = SwitchEndian(solbuf[6]);
    solheader_.iteration_ = SwitchEndian(solbuf[7]);
    solheader_.sumiteration_ = SwitchEndian(solbuf[8]);
    solheader_.numreact_ = SwitchEndian(solbuf[9]);
    solheader_.maxesz_ = SwitchEndian(solbuf[10]);
    solheader_.nummasters_ = SwitchEndian(solbuf[11]);
    solheader_.ptr_nodalsol_ = SwitchEndian(solbuf[12]);
    solheader_.ptr_elemsol_ = SwitchEndian(solbuf[13]);
    solheader_.ptr_react_ = SwitchEndian(solbuf[14]);
    solheader_.ptr_masters_ = SwitchEndian(solbuf[15]);
    solheader_.ptr_bc_ = SwitchEndian(solbuf[16]);
    solheader_.extrapolate_ = SwitchEndian(solbuf[17]);
    solheader_.mode_ = SwitchEndian(solbuf[18]);
    solheader_.symmetry_ = SwitchEndian(solbuf[19]);
    solheader_.complex_ = SwitchEndian(solbuf[20]);
    solheader_.numdofs_ = SwitchEndian(solbuf[21]);
    // jetzt Titel und Subtitel lesen
    memcpy(solheader_.title_, &solbuf[52], sizeof(int) * 20);
    solheader_.title_[79] = '\0';
    memcpy(solheader_.subtitle_, &solbuf[72], sizeof(int) * 20);
    solheader_.subtitle_[79] = '\0';
    // weiter gehts
    solheader_.changetime_ = SwitchEndian(solbuf[92]);
    solheader_.changedate_ = SwitchEndian(solbuf[93]);
    solheader_.changecount_ = SwitchEndian(solbuf[94]);
    solheader_.soltime_ = SwitchEndian(solbuf[95]);
    solheader_.soldate_ = SwitchEndian(solbuf[96]);
    solheader_.ptr_onodes_ = SwitchEndian(solbuf[97]);
    solheader_.ptr_oelements_ = SwitchEndian(solbuf[98]);
    solheader_.numexdofs_ = SwitchEndian(solbuf[99]);
    solheader_.ptr_extra_a_ = SwitchEndian(solbuf[100]);
    solheader_.ptr_extra_t_ = SwitchEndian(solbuf[101]);

    // DOFs einlesen
    delete[] solheader_.dof_;
    solheader_.dof_ = new int[solheader_.numdofs_];
    // erst mal die normalen DOFs kopieren
    for (i = 0; i < solheader_.numdofs_; ++i)
        solheader_.dof_[i] = SwitchEndian(solbuf[22 + i]);
    // exdofs reinhauen
    delete[] solheader_.exdof_;
    solheader_.exdof_ = new int[solheader_.numexdofs_];

    // Jetzt die TIME Varioable aus dem folgenden Daten besorgen
    // Zwei ints Ã¼berspringen
    fseek(rfp_, 2 * 4, SEEK_CUR);
    if (DoubleRecord(dsol, 100) != 100) {
        return (7);
    }
    // Time Variable sollte an erster Stelle stehen
    solheader_.time_ = dsol[0];

    offset = solheader_.offset_ + (solheader_.ptr_extra_a_ + 2) * 4;
    fseek(rfp_, (long)offset, SEEK_SET);
    int exdofbuf[64];
    if (IntRecord(exdofbuf, 64) != 64) {
        return (6);
    }
    for (i = 0; i < solheader_.numexdofs_; ++i)
        solheader_.exdof_[i] = exdofbuf[i];

    // buf_tab freed automatically when going out of scope
    return (0);
}

// -----------------------------------------------
// GetNodes
// -----------------------------------------------
// Liest die Koordinaten der Knoten ein
// Rueckgabe:
// 0        : alles ok
// 1        : Read Error Geometrieheader
// 2        : Read Error Nodes
// 3        : Read Error Elementbeschreibung
// 4        : Read Error ETYs
// 5        : Read Error Elementtabelle
// 6        : Read Error Elemente
int ReadRST::GetNodes(void)
{
    GeometryHeader ghdr;
    long int offset;
    int size, i, j;

    static const float DGR_TO_RAD = (float)M_PI / 180.0f;

    if (node_)
        return 0;

    // Springe erst mal zu der Geometrie-Tabelle
    offset = (rstheader_.ptr_geo_ + 1) * sizeof(int);
    fseek(rfp_, offset, SEEK_SET); // im pointer steht die Anzahl der int-Elemente vom Anfang

    size = 80;
    auto buf_up = std::make_unique<int[]>(size);

    if (fread(buf_up.get(), sizeof(int), size, rfp_) != size) {
        return (1);
    }

    // Werte zuweisen
    ghdr.maxety_ = SwitchEndian(buf_up[2]);
    ghdr.maxrl_ = SwitchEndian(buf_up[3]);
    ghdr.nodes_ = SwitchEndian(buf_up[4]);
    ghdr.elements_ = SwitchEndian(buf_up[5]);
    ghdr.maxcoord_ = SwitchEndian(buf_up[6]);
    ghdr.ptr_ety_ = SwitchEndian(buf_up[7]);
    ghdr.ptr_rel_ = SwitchEndian(buf_up[8]);
    ghdr.ptr_nodes_ = SwitchEndian(buf_up[9]);
    ghdr.ptr_sys_ = SwitchEndian(buf_up[10]);
    ghdr.ptr_eid_ = SwitchEndian(buf_up[11]);
    ghdr.ptr_mas_ = SwitchEndian(buf_up[16]);
    ghdr.coordsize_ = SwitchEndian(buf_up[17]);
    ghdr.elemsize_ = SwitchEndian(buf_up[18]);
    ghdr.etysize_ = SwitchEndian(buf_up[19]);
    ghdr.rcsize_ = SwitchEndian(buf_up[20]);

    // If a 32-bit field is zero, the 64-bit pointer parts may be stored
    // as hi/lo words elsewhere in buf_up. Combine hi/lo into a 64-bit value
    // depending on endianness handling.
    auto combine_pair = [&](int hi_idx, int lo_idx) -> long long {
        /*         if (SwitchEndian_ == DO_NOT_SWITCH) {
            return (static_cast<long long>(buf_up[hi_idx]) << 32) | static_cast<int>(buf_up[lo_idx]);
        } else {
            return (static_cast<long long>(SwitchEndian(buf_up[lo_idx])) << 32) |
                   static_cast<int>(SwitchEndian(buf_up[hi_idx]));
        } */
        if (SwitchEndian_ == DO_NOT_SWITCH) {
            return buf_up[hi_idx];
        } else {
            return SwitchEndian(buf_up[hi_idx]);
        }
    };

    if (ghdr.ptr_ety_ == 0) {
        std::cout << "64 bit first part: " << buf_up[21] << " second part: " << buf_up[22] << std::endl;
        ghdr.ptr_ety_ = combine_pair(21, 22);
    }
    if (ghdr.ptr_rel_ == 0) {
        std::cout << "64 bit first part: " << buf_up[23] << " second part: " << buf_up[24] << std::endl;
        ghdr.ptr_rel_ = combine_pair(23, 24);
    }
    if (ghdr.ptr_nodes_ == 0) {
        std::cout << "64 bit first part: " << buf_up[25] << " second part: " << buf_up[26] << std::endl;
        ghdr.ptr_nodes_ = combine_pair(25, 26);
    }
    if (ghdr.ptr_sys_ == 0) {
        std::cout << "64 bit first part: " << buf_up[27] << " second part: " << buf_up[28] << std::endl;
        ghdr.ptr_sys_ = combine_pair(27, 28);
    }
    if (ghdr.ptr_eid_ == 0) {
        std::cout << "64 bit first part: " << buf_up[29] << " second part: " << buf_up[30] << std::endl;
        ghdr.ptr_eid_ = combine_pair(29, 30);
    }
    if (ghdr.ptr_mas_ == 0) {
        std::cout << "64 bit first part: " << buf_up[31] << " second part: " << buf_up[32] << std::endl;
        ghdr.ptr_mas_ = combine_pair(31, 32);
    }

    std::cout << "Geometry Header Info:" << std::endl;
    std::cout << "first " << buf_up[0] << std::endl;
    std::cout << "unused position " << buf_up[1] << std::endl;
    std::cout << "Max ETYs: " << ghdr.maxety_ << std::endl;
    std::cout << "Max RELs: " << ghdr.maxrl_ << std::endl;
    std::cout << "Number of Nodes: " << ghdr.nodes_ << std::endl;
    std::cout << "Number of Elements: " << ghdr.elements_ << std::endl;
    std::cout << "Max Coordinates: " << ghdr.maxcoord_ << std::endl;
    std::cout << "Pointer to ETYs: " << ghdr.ptr_ety_ << std::endl;
    std::cout << "Pointer to RELs: " << ghdr.ptr_rel_ << std::endl;
    std::cout << "Pointer to Nodes: " << ghdr.ptr_nodes_ << std::endl;
    std::cout << "Pointer to System: " << ghdr.ptr_sys_ << std::endl;
    std::cout << "Pointer to EID: " << ghdr.ptr_eid_ << std::endl;
    std::cout << "Pointer to MAS: " << ghdr.ptr_mas_ << std::endl;
    std::cout << "Coordinate Size: " << ghdr.coordsize_ << std::endl;
    std::cout << "Element Size: " << ghdr.elemsize_ << std::endl;
    std::cout << "ETY Size: " << ghdr.etysize_ << std::endl;
    std::cout << "RC Size: " << ghdr.rcsize_ << std::endl;

    // Jetzt zu den KNoten springen und diese lesen (Lead in ueberspringen)
    offset = (ghdr.ptr_nodes_ + PTR_OFFSET) * sizeof(int);
    fseek(rfp_, (long)offset, SEEK_SET);

    // Jetzt die NODES definieren
    node_ = new Node[ghdr.nodes_];
    size_t NODE_SIZE = 7;
    for (i = 0; i < ghdr.nodes_; ++i) {
        double nodeInfo[NODE_SIZE];
        if (DoubleRecord(nodeInfo, NODE_SIZE) != NODE_SIZE)
            return (2);

        // Werte jetzt umdrehen
        node_[i].id_ = nodeInfo[0];
        node_[i].x_ = nodeInfo[1];
        node_[i].y_ = nodeInfo[2];
        node_[i].z_ = nodeInfo[3];
        node_[i].thxy_ = DGR_TO_RAD * nodeInfo[4];
        node_[i].thyz_ = DGR_TO_RAD * nodeInfo[5];
        node_[i].thzx_ = DGR_TO_RAD * nodeInfo[6];
        node_[i].MakeRotation();
        if (mode64_) {
            fseek(rfp_, 12, SEEK_CUR);
        }
    }
    std::cout << "First 10 Nodes:" << std::endl;
    for (j = 0; j < 10 && j < ghdr.nodes_; ++j) {
        std::cout << "Node " << node_[j].id_ << ": " << node_[j].x_ << ", " << node_[j].y_ << ", " << node_[j].z_
                  << std::endl;
        std::cout << "  Rotation angles (rad): " << node_[j].thxy_ << ", " << node_[j].thyz_ << ", " << node_[j].thzx_
                  << std::endl;
    }
    // das war's!
    anznodes_ = ghdr.nodes_;

    // Jetzt die Elemente lesen: zuerst mal zu ETY um die Elementbeschreibungen zu laden
    offset = (ghdr.ptr_ety_ + 2) * 4;
    fseek(rfp_, (long)offset, SEEK_SET);

    // buf_up still in scope, reuse if appropriate or create new buffer:
    auto buf2 = std::make_unique<int[]>(ghdr.maxety_);
    if (fread(buf2.get(), sizeof(int), ghdr.maxety_, rfp_) != ghdr.maxety_) {
        return (3);
    }

    // ETYs im Objekt erstellen
    ety_.resize(ghdr.maxety_);
    anzety_ = ghdr.maxety_;
    auto etybuf_up = std::make_unique<int[]>(ghdr.etysize_);

    for (i = 0; i < ghdr.maxety_; ++i) {
        if (mode64_) {
            offset = (SwitchEndian(buf2[i]) + ghdr.ptr_ety_ + 2) * 4; //TODO: is offset correct?
        } else {
            offset = (SwitchEndian(buf2[i]) + 2) * 4;
        }
#ifdef WIN32
        _fseeki64(rfp_, offset, SEEK_SET);
#else
        fseek(rfp_, offset, SEEK_SET);
#endif

        if (fread(etybuf_up.get(), sizeof(int), ghdr.etysize_, rfp_) != ghdr.etysize_) {
            return (4);
        }
        // Daten jetzt in den ety bringen
        ety_[i].id_ = SwitchEndian(etybuf_up[0]);
        ety_[i].routine_ = SwitchEndian(etybuf_up[1]);
        ety_[i].keyops_[0] = SwitchEndian(etybuf_up[2]);
        ety_[i].keyops_[1] = SwitchEndian(etybuf_up[3]);
        ety_[i].keyops_[2] = SwitchEndian(etybuf_up[4]);
        ety_[i].keyops_[3] = SwitchEndian(etybuf_up[5]);
        ety_[i].keyops_[4] = SwitchEndian(etybuf_up[6]);
        ety_[i].keyops_[5] = SwitchEndian(etybuf_up[7]);
        ety_[i].keyops_[6] = SwitchEndian(etybuf_up[8]);
        ety_[i].keyops_[7] = SwitchEndian(etybuf_up[9]);
        ety_[i].keyops_[8] = SwitchEndian(etybuf_up[10]);
        ety_[i].keyops_[9] = SwitchEndian(etybuf_up[11]);
        ety_[i].keyops_[10] = SwitchEndian(
            etybuf_up[12]); //TODO: is most negative number possible in file.rst: does this make sense? offset to big?
        ety_[i].keyops_[11] = SwitchEndian(etybuf_up[13]); //28
        ety_[i].dofpernode_ = SwitchEndian(etybuf_up[33]); //89
        ety_[i].nodes_ = SwitchEndian(etybuf_up[60]);
        ety_[i].nodeforce_ = SwitchEndian(etybuf_up[62]);
        ety_[i].nodestress_ = SwitchEndian(etybuf_up[93]);
    }
    // Passt schon!
    // Jetzt die Elemente selber einlesen
    element_.resize(ghdr.elements_);
    anzelem_ = ghdr.elements_;

    // hinsurfen und lesen
    offset = (ghdr.ptr_eid_ + 2) * 4;
#ifdef WIN32
    _fseeki64(rfp_, offset, SEEK_SET);
#else
    fseek(rfp_, offset, SEEK_SET);
#endif

    // allocate buf3 for element table
    auto buf3 = std::make_unique<int[]>(2 * ghdr.elements_);
    if (fread(buf3.get(), sizeof *buf3.get(), 2 * ghdr.elements_, rfp_) != 2 * ghdr.elements_) {
        std::cerr << "Problems occurred while reading the element table\n";
        return 5;
    }

    auto etybuf2 = std::make_unique<int[]>(10);
    // Element accounting for user information
    int populations[ANSYS::LIB_SIZE];
    memset(populations, 0, sizeof(int) * ANSYS::LIB_SIZE);
    // Jetzt mit Schleife alle Elemente packen
    for (i = 0; i < ghdr.elements_; ++i) {
        if (mode64_) {
            offset = (SwitchEndian(static_cast<uint32_t>(buf3[2 * i])) + ghdr.ptr_eid_ + 2) * 4;
        } else {
            offset = (SwitchEndian(static_cast<uint32_t>(buf3[i])) + 2) * 4;
        }

#ifdef WIN32
        _fseeki64(rfp_, offset, SEEK_SET);
#else
        fseek(rfp_, offset, SEEK_SET);
#endif
        if (fread(etybuf2.get(), sizeof(int), 10, rfp_) != 10) { //TODO: invalid first index in type_
            return (6);
        }
        // Jetzt Daten zuweisen
        element_[i].material_ = SwitchEndian(etybuf2[0]);
        element_[i].type_ = SwitchEndian(etybuf2[1]);
        element_[i].real_ = SwitchEndian(etybuf2[2]);
        element_[i].section_ = SwitchEndian(etybuf2[3]);
        element_[i].coord_ = SwitchEndian(etybuf2[4]);
        element_[i].death_ = SwitchEndian(etybuf2[5]);
        element_[i].solidmodel_ = SwitchEndian(etybuf2[6]);
        element_[i].shape_ = SwitchEndian(etybuf2[7]);
        element_[i].num_ = SwitchEndian(etybuf2[8]);
        element_[i].anznodes_ = 0;

        for (j = 0; j < anzety_; ++j) {
            if (ety_[j].id_ == element_[i].type_) // find ety position -> make faster
            {
                element_[i].anznodes_ = ety_[j].nodes_;
                ++(populations[ety_[j].routine_]);
            }
        }
        element_[i].nodes_ = new int[element_[i].anznodes_];
        if (IntRecord(element_[i].nodes_, element_[i].anznodes_) != element_[i].anznodes_) {
            return (6);
        }
    }

    // Some statistics
    std::cout << "Total number of nodes: " << ghdr.nodes_ << std::endl;
    std::cout << "Total number of elements: " << ghdr.elements_ << std::endl;

    // std::string InfoMessage("Number of elements in each element category:");
    for (i = 0; i < ANSYS::LIB_SIZE; ++i) {
        if (populations[i]) {
            // char numBuf[128];
            // sprintf(numBuf," in routine %d, %d elements;",i,populations[i]);
            // InfoMessage += numBuf;

            std::cout << "Element type [" << i << "]: " << populations[i] << " elements" << std::endl;
        }
    }
    // size_t pos = InfoMessage.find_last_of(";");
    // if(pos != std::string::npos)
    //   InfoMessage.resize(pos-1);
    // InfoMessage += ".";
    // sendInfo("%s", InfoMessage.c_str());

    // buf3 and etybuf2 freed automatically
    return (0);
}

// -----------------------------------------------
// GetTime
// -----------------------------------------------
// Return:
// <0   Fehler und zwar:
// -1.0     : Index out of Range
// -2.0     : Keine Zeittafel vorhanden
double ReadRST::GetTime(int pos)
{
    if (pos > rstheader_.numsets_) // index ausserhalb
        return (-1.0);

    if (timetable_ == NULL) // keine Zeittabelle gelesen
    {
        return pos;
        // return(-2.0);
    }

    return (timetable_[pos]); // Zeitwert zurÃ¼ck
}

// -----------------------------------------------
// ReadFortranRecordInts
// -----------------------------------------------
// helper function to read Fortran-style records with lead/trailer integers
int ReadRST::ReadFortranRecordInts(int32_t *buf, size_t expectedCount, int headWords)
{
    if (!rfp_)
        return -1;

    if (headWords > 0) {
        std::vector<int32_t> lead(headWords);
        if (fread(lead.data(), sizeof(int32_t), headWords, rfp_) != (size_t)headWords)
            return -2;
    }
    if (expectedCount) {
        if (buf) {
            if (fread(buf, sizeof(int32_t), expectedCount, rfp_) != expectedCount)
                return -3;
        } else {
            // no buffer provided: skip payload
            if (fseek(rfp_, static_cast<long>(expectedCount * sizeof(int32_t)), SEEK_CUR) != 0)
                return -3;
        }
    }
    /* 
    std::vector<int32_t> trail(headWords);
    if (fread(trail.data(), sizeof(int32_t), headWords, rfp_) != (size_t)headWords)
        return -4;

    if (lead != trail) {
        std::clog << "Warning: Fortran record lead/trail mismatch (lead != trail)\n";
    } */

    return static_cast<int>(expectedCount);
}
