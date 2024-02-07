/* XDMF_out.h */
/* Wolfgang Tichy 2/2024 */

/* VisIt and ParaView can read XDMF format.
   Our output in XDMF format consists of a .xml file and a .bin file. The
   actual data is in the .bin file, while the .xml file only contains a
   description of the data in XML format.

   Some description about XDMF is on
   http://www.xdmf.org
   http://www.xdmf.org/index.php/XDMF_Model_and_Format
   http://www.paraview.org/Wiki/ParaView/Data_formats

   What we have here is quite similar to what is in bamps 4.0. But I added
   "Time Value" and "Spatial" to B_spatial. Now VisIt shows the correct time
   and paraview can show more than t=0. */


/* XML format strings to make .xmf files using fprintf,
   based on bamps and https://www.paraview.org/Wiki/ParaView/Data_formats */
static const char *B_head_xmf =
  "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
  "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.1\">\n"
  "  <Domain>\n";
static const char *E_head_xmf =
  "  </Domain>\n"
  "</Xdmf>\n";

static const char *B_temporal_xmf =
  "    <Grid CollectionType=\"Temporal\" GridType=\"Collection\" Name=\"TCollection\">\n"
  "      <Geometry Type=\"None\"/>\n"
  "      <Topology Dimensions=\"0\" Type=\"NoTopology\"/>\n";
static const char *E_temporal_xmf =
  "    </Grid>\n";

static const char *B_spatial_xmf =
  "      <Grid CollectionType=\"Spatial\" GridType=\"Collection\" Name=\"SCollection\">\n"
  "        <Time Value=\"%.9lf\"/>\n"
  "        <Geometry Type=\"None\"/>\n"
  "        <Topology Dimensions=\"0\" Type=\"NoTopology\"/>\n";
static const char *E_spatial_xmf =
  "      </Grid>\n";

static const char *B_E_grid_xmf =
  "        <Grid Name=\"%s\">\n"
  "          <Time Value=\"%.9lf\"/>\n"
  "          <Geometry Type=\"XYZ\">\n"
  "            <DataItem DataType=\"Float\" Dimensions=\"%d %d\" Format=\"%s\" Seek=\"%ld\" Precision=\"4\">\n"
  "              %s\n"
  "            </DataItem>\n"
  "          </Geometry>\n"
  "          <Topology Dimensions=\"%d %d %d\" Type=\"3DSMesh\"/>\n"
  "          <Attribute Center=\"Node\" Name=\"%s\" Type=\"Scalar\">\n"
  "            <DataItem DataType=\"Float\" Dimensions=\"%d %d %d\" Format=\"%s\" Seek=\"%ld\" Precision=\"4\">\n"
  "              %s\n"
  "            </DataItem>\n"
  "          </Attribute>\n"
  "        </Grid>\n";

/* same format strings, but in WT's simplified text format */
static const char *B_head_smf =
  "# sdmf:   binarydata: float   TopologyType: 3DSMesh   AttributeCenter: Node\n"
  "# node n[0] n[1] n[2] xyzseek varseek\n";
static const char *E_head_smf = "";
static const char *B_temporal_smf = "";
static const char *E_temporal_smf = "";
static const char *B_spatial_smf = "\n# \"time = %.16lg\"\n";
static const char *E_spatial_smf = "";
static const char *B_E_grid_smf = "%s\t%d %d %d\t%ld %ld\n";
