#ifndef REACTION_H
#define REACTION_H


#include <math.h>
#include <iostream>
#include <string>
using namespace std;


class reaction
{
 private:
  
  const static int maxNoReactionSpecies = 4;

  double reactionRate;

  int reactantSpeciesList[maxNoReactionSpecies+1];
  int productSpeciesList[maxNoReactionSpecies+1];
  
  double reactionRateFunction(int j, double Tgas, double Telectron);
  
  
 public:
  
  reaction(void);

  void setReactionRate(int j, double Tgas, double Telectron);
  void setReactantAndProductSpecies(int j);
  
  int returnNumberOfReactants(void);
  int returnNumberOfProducts(void);
  
  int returnReactant(int);
  int returnProduct(int);
  
  double returnReactionRate();
  
};
#endif
