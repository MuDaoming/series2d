#ifndef SETDES_HPP
#define SETDES_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <map>

#include "Intg.hpp"
#include "linear.hpp"


// // --------------------------------Branch　Functions-----------------------------

bool fullbranch(std::vector<int>& info, std::vector<int>& index, int B) {
  std::vector<int> sum(B, 0);

  for (int i = 0; i < index.size(); i++) {
    if (index[i] == 1) sum[info[i]] += 1;
  }

  for (int b = 0; b < B; b++) {
    if (sum[b] == 0) return false;
  }

  return true;
}

template <typename FT>
std::vector<int> branchinfo(std::vector<std::vector<FT>>& topS, int B) {
  int len = topS.size();
  std::vector<int> info(len - B);

  for (int i = 0; i < B; i++) {
    for (int j = B; j < len; j++) {
      if (isMagnitudeBelowEpsilon(topS[i][j] - (FT)1) == 1)
        info[j - B] = i;  // caution here!
    }
  }

  return info;
}

// // --------------------------------- MIs Functions ---------------------------------

template <typename FT>
void sortMis(std::vector<Intg<FT>>& mis) {
  int numOfMis = mis.size();
  for (int i = 0; i < numOfMis; i++) {
    for (int j = i + 1; j < numOfMis; j++) {
      if (mis[j] < mis[i]) {
        std::swap(mis[i], mis[j]);
      }
    }
  }
}

template <typename FT>
int findPositionOfMi(Intg<FT>& mi, std::vector<Intg<FT>>& mis, int num) {
  for (int i = num - 1; i >= 0; i--) {
    if (mi == mis[i]) {
      return i;
    }
  }
  return -1;
}

template <typename FT>
void dropPropagatorInMi(Intg<FT>& mi, int i,int L)
// will change input from mi to submi
// input should be a corner integral
{
  int count = -1;
  for (int j = 0; j < mi.index.size(); j++) {
    if (mi.index[j] > 1 || mi.index[j] < 0) {
      // mi is not a corner integral
      throw std::invalid_argument("Error: mi is not a corner integral.");
    }
    if (mi.index[j] == 1) {
      count++;
    }
    if (count == i) {
      mi.index[j] = mi.index[j] - 1;
      mi.d = mi.d - (FT)2/(FT)L;
      break;
    }
  }
}

// // --------------------------------- Set DEs Functions ---------------------------------

template <typename FT>
std::vector<std::vector<FT>> getSubS(std::vector<std::vector<FT>>& topS,
                                     std::vector<int>& index, int B) {
  int n = topS.size();
  int nu = index.size();

  // dimension of sub S
  int sub_n = B;
  for (int i = 0; i < nu; i++) {
    if (index[i] == 1) {
      sub_n++;
    }
  }

  // create sub S
  std::vector<std::vector<FT>> subS(sub_n, std::vector<FT>(sub_n));

  int subrow = 0;
  for (int i = 0; i < n; i++) {
    if (i >= B && index[i - B] != 1) continue;  // only if the index equals 1
    int subcol = 0;
    for (int j = 0; j < n; j++) {
      if (j >= B && index[j - B] != 1) continue;  // only if the index equals 1
      subS[subrow][subcol] = topS[i][j];
      subcol++;
    }
    subrow++;
  }

  return subS;
}

int binary(const std::vector<int>& index) {
  int total = 0;
  int len = index.size();
  for (int i = 0; i < len; i++) total += index[i] * pow(2, len - i - 1);
  return total;
}



// check there is no pole in the fourth quadrant
template <typename FT>
int checkAnyPoleNotInFourthQuadrant(std::vector<FT>& poles, std::vector<int>& IsPoleHere) {
  // first check if type FT is complex-like using is_complex_like_v
    if constexpr (is_complex_like_v<FT>) {
      for (int i = 0; i < poles.size(); i++) {
        if (IsPoleHere[i] == 1 && poles[i].imag() < 0 && poles[i].real() > 0) {
          throw std::runtime_error("Error: There is a pole in the fourth quadrant, which is not allowed.");
          return 0; // if any pole is in the fourth quadrant, return 0
        }
      }
      return 1; // if all poles are not in the fourth quadrant, return 1
    } else {
      return 1; // if all poles are real
    } 
}


template <typename FT>
void setDes(std::vector<std::vector<FT>>& topS, int B, int P, int L, FT d0,
            std::vector<Intg<FT>>& mis, std::vector<int>& IsPoleHere,
            std::vector<FT>& poles, std::vector<std::vector<FT>>& coematrix) {

  // find MIs
  int total = 1 << P;
  Intg<FT> temp_mis[total];
  int numofmis = 0;
  std::map<int, std::vector<FT>> sol;

  std::vector<int> info = branchinfo(topS, B);

  for (int i = total - 1; i > 0; i--) {
    std::vector<int> index(P);
    int num = i;
    int nu = 0;
    // Fill the index based on its binary representation
    for (int j = P - 1; j >= 0; j--) {
      index[j] = num % 2;
      nu += index[j];
      num /= 2;
    }

    if (!fullbranch(info, index, B)) continue;

    std::vector<std::vector<FT>> S = getSubS(topS, index, B);

    if (IsFullRank(S) == 1) {
      std::vector<FT> z(P + B, (FT)0);  // modified to first B element
      for (int b = 0; b < B; b++) z[b] = (FT)1;
      solve_linear(S, z);
      sol[i] = z;
    } else {
      std::vector<FT> z(P + B, (FT)0);
      findNullSpace(S, z, B);
      if (isMagnitudeBelowEpsilon(sumB(z, B)) == 1) {  // caution here!
        continue;
      } else {
        sol[i] = z;
      }
    }

    numofmis++;
    Intg<FT> mi = {d0 - (FT)(2 * (P - nu)) / (FT)L, index};
    temp_mis[numofmis - 1] = mi;
  }

  mis.clear();
  for (int i = 0; i < numofmis; i++) mis.push_back(temp_mis[i]);

  if (numofmis == 0) {
    throw std::invalid_argument(
        "Error: no MI found under the given kinematics, all the integrals have "
        "vanishing det(R) and det(S).");
  }

  // sort
  sortMis(mis);

  IsPoleHere.resize(numofmis);
  poles.resize(numofmis);
  coematrix.resize(numofmis, std::vector<FT>(numofmis));

  for (int i = 0; i < numofmis; i++) {
    Intg<FT> mi = mis[i];

    int nu = 0;
    for (int in = 0; in < mi.index.size(); in++) nu += mi.index[in];

    std::vector<std::vector<FT>> S = getSubS(topS, mi.index, B);

    if (nu == -1) {  // Need to be deleted
      IsPoleHere[i] = 1;
      poles[i] = - S[1][1] / (FT)2;
      coematrix[i][i] = (mi.d - (FT)1 - (FT)nu) / (FT)2;
    } else {
      if (1 == IsFullRank(S)) {
        IsPoleHere[i] = 1;

        std::vector<FT> z = sol[binary(mi.index)];

        // if (1 == isMagnitudeBelowEpsilon(sumB(z, B))) 
        // {
        //   poles[i] = (FT)0;
        // } else {
        //   poles[i] = sumB(z, B) / (FT)2;  // caution here!
        // }

        poles[i] = sumB(z, B) / (FT)2;   // caution here!!!

        coematrix[i][i] = ((FT)L * mi.d - (FT)B - (FT)nu) / (FT)2;

        for (int j = 0; j < nu; j++) {
          // if (1 == isMagnitudeBelowEpsilon(z[j + B])) {
          //   continue;
          // }
          Intg<FT> submi = mi;
          dropPropagatorInMi(submi, j,L);
          if (fullbranch(info, submi.index, B)) {
            int pos = findPositionOfMi(submi, mis, numofmis);
            // if (pos == -1) {
            //   throw std::runtime_error(
            //       "Error: master integrals are incomplete, there is a integral "
            //       "with vanishing det(S) and det(R) in differential "
            //       "equations.");
            // }
            if (pos == -1) {
              if(!isMagnitudeBelowEpsilon(z[j + B]) ) {
              std::cerr << "Warning: a z-coefficient supposed to be zero is " << std::setprecision(3) <<z[j + B] << std::endl;}
              continue;
              }
            coematrix[i][pos] = z[j + B] / (FT)2;
          }
        }
      } else {
        IsPoleHere[i] = 0;

        std::vector<FT> z = sol[binary(mi.index)];

        for (int j = 0; j < nu; j++) {
          // if (isMagnitudeBelowEpsilon(z[j + B])) {
          //   continue;
          // }
          Intg<FT> submi = mi;
          dropPropagatorInMi(submi, j,L);
          if (fullbranch(info, submi.index, B)) {
            int pos = findPositionOfMi(submi, mis, numofmis);
            // if (pos == -1) {
            //   throw std::runtime_error(
            //       "Error: master integrals are incomplete, there is a integral "
            //       "with vanishing det(S) and det(R) in differential "
            //       "equations.");
            // }
            if (pos == -1) {
              if(!isMagnitudeBelowEpsilon(z[j + B]) ) {
              std::cerr << "Warning: a z-coefficient supposed to be zero is " << std::setprecision(3) <<z[j + B] << std::endl;}
              continue;
              }
            coematrix[i][pos] = -z[j + B] / sumB(z, B);
          }
        }
      }
    }
  }

  // check all the poles are not in the fourth quadrant
  checkAnyPoleNotInFourthQuadrant(poles, IsPoleHere);
}



// -----------------------------Functions After Setting DEs------------------------------------------


template <typename FT>
std::vector<FT> getFarestPoles(
    std::vector<int>& IsPoleHere, std::vector<FT>& poles,
    std::vector<std::vector<FT>>&
        coematrix)  // find farest poles for sorted DEs
{
  int numOfMis = IsPoleHere.size();
  std::vector<FT> farestpoles(numOfMis);
  for (int i = 0; i < numOfMis; i++) {
    farestpoles[i] = IsPoleHere[i] ? absSafe(poles[i]) : (FT)0;
    for (int j = 0; j < i; j++) {
      if (!isMagnitudeBelowEpsilon(coematrix[i][j]) &&
          absSafe(farestpoles[j]) > absSafe(farestpoles[i])) {
        farestpoles[i] = farestpoles[j];
      }
    }
  }
  return farestpoles;
}

template <typename FT, typename CT>
void getBcs(int B, std::vector<Intg<FT>>& mis, std::vector<CT>& bc, std::vector<std::vector<FT>>& topS , int L) {
  std::vector<int> info = branchinfo(topS, B);
  int numOfMis = mis.size();
  int n = ((mis[0]).index).size();
  bc.resize(numOfMis);

  int nu0 = 0;
  for (int i = 0; i < n; i++) nu0 += (mis[0]).index[i];
  
  FT num_gamma = boost::math::tgamma(real((CT)nu0 - (CT)L*(CT)mis[0].d / (CT)2));
  // note that we use the condition that all the integrals in DEs have the same
  // regions at infinity



  // if we do not use numOfMis but mis.size() directly, under Ofast,
  // here would be a segmentation fault
  for (int i = 0; i < numOfMis; i++) {

    std::vector<int> sum(B, 0);
    for (int j = 0; j < n; j++) {
      if (mis[i].index[j] == 1) sum[info[j]] += 1;
    }

    bc[i] = num_gamma;
    for (int b = 0; b < B; b++) bc[i] /= boost::math::tgamma(sum[b])*pow(-1,sum[b]);
  }
}

template <typename FT, typename CT>
std::vector<int> getPositionsOfInfCutsInRegPoints(std::vector<FT>& farestpoles,
                                                  std::vector<CT>& run) {
  std::vector<int> position_of_inf_cut(farestpoles.size());
  for (int i = 0; i < farestpoles.size(); i++) {
    for (int j = run.size() - 1; j >= 0; j--) {
      if (absSafe(run[j]) > 1.95 * absSafe(farestpoles[i])) {
        position_of_inf_cut[i] = j;
        break;
      }
    }
  }
  return position_of_inf_cut;
}

#endif