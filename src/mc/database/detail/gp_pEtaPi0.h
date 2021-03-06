#include "ProductionDataBase.h"

using namespace ant;
using namespace ant::mc::data;


static ProductionDataBase::XSections_t::value_type gp_pEtaPi0 =
{ ParticleTypeTreeDatabase::Channel::gp_pEtaPi0,
  ProductionDataBase::MakeInterPolator({
    
          { 931,   0.0},  // threshold
          // cb/cbElsa data
{937,   0.0010 },
{949,   0.0020 },
{961,   0.0060 },
{972,   0.012 },
{983,   0.025 },
{995,   0.045 },
{1006,   0.071 },
{1017,   0.107 },
{1028,   0.144 },
{1039,   0.192 },
{1050,   0.243 },
{1061,   0.33 },
{1072,   0.401 },
{1082,   0.491 },
{1093,   0.573 },
{1103,   0.664 },
{1114,   0.773 },
{1124,   0.885 },
{1134,   1.061 },
{1144,   1.149 },
{1154,   1.24 },
{1164,   1.473 },
{1174,   1.548 },
{1184,   1.594 },
{1193,   1.731 },
{1203,   1.904 },
{1212,   2.105 },
{1221,   2.076 },
{1230,   2.158 },
{1239,   2.265 },
{1248,   2.289 },
{1257,   2.406 },
{1266,   2.519 },
{1275,   2.659 },
{1283,   2.718 },
{1291,   2.666 },
{1300,   2.683 },
{1308,   2.904 },
{1316,   3.035 },
{1322,   2.941 },
{1335,   3.127 },
{1343,   3.025 },
// fit to cb/cbElsa data: f(x) = p0 + p1 * x + p1 * x^2
//  Minimizer is Linear
//  Chi2                      =      1.62506
//  NDf                       =           17
//  p0                        =     -33.1171   +/-   18.2706     
//  p1                        =      46.9949   +/-   25.0073     
//  p2                        =     -14.8619   +/-   8.5167  
{1350,   3.24022 },
{1356,   3.28089 },
{1362,   3.3205 },
{1368,   3.35903 },
{1374,   3.39649 },
{1380,   3.43288 },
{1386,   3.4682 },
{1392,   3.50246 },
{1398,   3.53564 },
{1404,   3.56775 },
{1410,   3.59879 },
{1416,   3.62876 },
{1422,   3.65766 },
{1428,   3.68549 },
{1434,   3.71226 },
{1440,   3.73795 },
{1446,   3.76257 },
{1452,   3.78612 },
{1458,   3.8086 },
{1464,   3.83001 },
{1470,   3.85035 },
{1476,   3.86962 },
{1482,   3.88783 },
{1488,   3.90496 },
{1494,   3.92102 },
{1500,   3.93601 },
{1506,   3.94993 },
{1512,   3.96278 },
{1518,   3.97456 },
{1524,   3.98527 },
{1530,   3.99491 },
{1536,   4.00348 },
{1542,   4.01098 },
{1548,   4.01741 },
{1554,   4.02277 },
{1560,   4.02706 },
{1566,   4.03028 },
{1572,   4.03243 },
{1578,   4.03351 },
{1584,   4.03352 },
{1590,   4.03246 },
{1596,   4.03033 },
{1602,   4.02713 },
{1608,   4.02286 },
{1614,   4.01752 },
{1620,   4.01111 },
{1626,   4.00363 },
{1632,   3.99508 },
{1638,   3.98546 },
{1644,   3.97477 }

  })};

