#include <dolfin/elements/FE.h>

namespace FE
{

//-----------------------------------------------------------------------------
Array<std::string> const init_families()
{
  Array<std::string> ret;
  ret.push_back(LAGRANGE);
  ret.push_back(DG);
  ret.push_back(BDM);
  return ret;
}

//-----------------------------------------------------------------------------
Array<std::string> const init_elements()
{
  Array<std::string> ret;
  ret.push_back(LAGRANGE1DP1S);
  ret.push_back(LAGRANGE1DP2S);
  ret.push_back(LAGRANGE2DP1S);
  ret.push_back(LAGRANGE2DP2S);
  ret.push_back(LAGRANGE3DP1S);
  ret.push_back(LAGRANGE3DP2S);
  ret.push_back(LAGRANGE2DP1V);
  ret.push_back(LAGRANGE2DP2V);
  ret.push_back(LAGRANGE3DP1V);
  ret.push_back(LAGRANGE3DP2V);
  ret.push_back(DG1DP0S);
  ret.push_back(DG1DP1S);
  ret.push_back(DG1DP2S);
  ret.push_back(DG2DP0S);
  ret.push_back(DG2DP1S);
  ret.push_back(DG2DP2S);
  ret.push_back(DG3DP0S);
  ret.push_back(DG3DP1S);
  ret.push_back(DG3DP2S);
  ret.push_back(DG2DP0V);
  ret.push_back(DG2DP1V);
  ret.push_back(DG2DP2V);
  ret.push_back(DG3DP0V);
  ret.push_back(DG3DP1V);
  ret.push_back(DG3DP2V);
  ret.push_back(BDM1DP2);
  return ret;
}

}
