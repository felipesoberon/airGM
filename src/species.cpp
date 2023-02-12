#include "species.h"


species::species(void)
{
  int i;
  for (i=0;i<=maxNoReactions;i++)
    {
      listOfReactionNoSources[i] = 0;
      listOfReactionNoLosses[i] = 0;
      listOfReactionNoSourcesMultiplier[i] = 0;
      listOfReactionNoLossesMultiplier[i] = 0;
    }
}


void species::setFormula(int speciesIndex)
{
  switch(speciesIndex)
  {
    case 0:
      formula="M";
      break;
    case 1:
      formula="N+";
      break;
    case 2:
      formula="N2+";
      break;
    case 3:
      formula="N3+";
      break;
    case 4:
      formula="N4+";
      break;
    case 5:
      formula="O+";
      break;
    case 6:
      formula="O2+";
      break;
    case 7:
      formula="O4+";
      break;
    case 8:
      formula="NO+";
      break;
    case 9:
      formula="N2O+";
      break;
    case 10:
      formula="NO2+";
      break;
    case 11:
      formula="H+";
      break;
    case 12:
      formula="H2+";
      break;
    case 13:
      formula="H3+";
      break;
    case 14:
      formula="OH+";
      break;
    case 15:
      formula="H2O+";
      break;
    case 16:
      formula="H3O+";
      break;
    case 17:
      formula="e";
      break;
    case 18:
      formula="O-";
      break;
    case 19:
      formula="O2-";
      break;
    case 20:
      formula="O3-";	
      break;
    case 21:
      formula="O4-";	
      break;
    case 22:
      formula="NO-";	
      break;
    case 23:
      formula="N2O-";	
      break;
    case 24:
      formula="NO2-";	
      break;
    case 25:
      formula="NO3-";	
      break;
    case 26:
      formula="H-";break;
    case 27:
      formula="OH-";break;
    case 28:
      formula="N(2_D)";break;
    case 29:
      formula="N2(A_3_Sigma)";break;
    case 30:
      formula="N2(B_3_Pi)";break;
    case 31:
      formula="O(1_D)";break;
    case 32:
      formula="H";break;
    case 33:
      formula="N";break;
    case 34:
      formula="O";break;
    case 35:
      formula="O2(a_1_Delta)";break;
    case 36:
      formula="O3";break;
    case 37:
      formula="NO";break;
    case 38:
      formula="N2O";break;
    case 39:
      formula="NO2";break;
    case 40:
      formula="NO3";break;
    case 41:
      formula="N2O3";break;
    case 42:
      formula="N2O4";break;
    case 43:
      formula="N2O5";break;
    case 44:
      formula="H2";break;
    case 45:
      formula="OH";break;
    case 46:
      formula="HO2";break;
    case 47:
      formula="H2O2";break;
    case 48:
      formula="HNO";break;
    case 49:
      formula="HNO2";break;
    case 50:
      formula="HNO3";break;
    case 51:
      formula="N2";break;
    case 52:
      formula="O2";break;
    case 53:
      formula="H2O";break;
    default:
      formula="X";
    }
}


void species::setDensity(double densityvalue)
{
  if (densityvalue < minimumDensity)
    density = 0;
  else
    density = densityvalue;
}



void species::setLoss(double lossvalue)
{ loss = lossvalue; }

void species::setSource(double sourcevalue)
{ source = sourcevalue; }



void species::incrementNumberOfReactionNoSources(void)
{ listOfReactionNoSources[0]++; }

void species::incrementNumberOfReactionNoLosses(void)
{ listOfReactionNoLosses[0]++; } 

void species::setReactionNoSourcesItem(int index, int value)
{ listOfReactionNoSources[index] = value; }

void species::setReactionNoLossesItem(int index, int value)
{ listOfReactionNoLosses[index] = value; }

void species::setReactionNoSourcesMultiplierItem(int index, int value)
{ listOfReactionNoSourcesMultiplier[index] = value; }

void species::setReactionNoLossesMultiplierItem(int index, int value)
{ listOfReactionNoLossesMultiplier[index] = value; }



int species::returnNumberOfReactionNoSources(void)
{ return listOfReactionNoSources[0]; }

int species::returnNumberOfReactionNoLosses(void)
{ return listOfReactionNoLosses[0]; }


int species::returnReactionNoSources(int index)
{ return listOfReactionNoSources[index]; }

int species::returnReactionNoLosses(int index)
{ return listOfReactionNoLosses[index]; }


int species::returnReactionNoSourcesMultiplier(int index)
{ return listOfReactionNoSourcesMultiplier[index]; }

int species::returnReactionNoLossesMultiplier(int index)
{ return listOfReactionNoLossesMultiplier[index]; }



double species::returnDensity(void)
{ return density; }

double species::returnLoss(void)
{ return loss; }

double species::returnSource(void)
{ return source; }

string species::returnFormula(void)
{ return formula; }


void species::printListOfSourcesAndLossesReactions(void)
{
  int i;
  int numberOfSources = listOfReactionNoSources[0];
  int numberOfLosses  = listOfReactionNoLosses[0];
  int reactionNumber;
  int reactionMultiplier;

  cout << formula << ",sources:,";
  cout << numberOfSources;
  cout << ",";

  for (i=1; i<=numberOfSources; i++)
    {
      reactionNumber     = listOfReactionNoSources[i];
      reactionMultiplier = listOfReactionNoSourcesMultiplier[i];
      cout << reactionNumber;
      cout << "(x";
      cout << reactionMultiplier;
      cout << "),";
    }
  cout << endl;


  cout << formula << ",losses:,";
  cout << numberOfLosses;
  cout << ",";

  for (i=1; i<=numberOfLosses; i++)
    {
      reactionNumber     = listOfReactionNoLosses[i];
      reactionMultiplier = listOfReactionNoLossesMultiplier[i];
      cout << reactionNumber;
      cout << "(x";
      cout << reactionMultiplier;
      cout << "),";
    }
  cout << endl;
  
}







void species::processReduceListsOfReactionNo(void)
{
  int is, il;
  int NoReactionsSources = listOfReactionNoSources[0];
  int NoReactionsLosses = listOfReactionNoLosses[0];
  int lossReactionNo;
  int sourceReactionNo;
  int lossReactionNoMultiplier;
  int sourceReactionNoMultiplier;

  for (is=1; is<=NoReactionsSources; is++)
    {
      sourceReactionNo = listOfReactionNoSources[is];
      
      for (il=1; il<=NoReactionsLosses; il++)
	{
	  lossReactionNo = listOfReactionNoLosses[il];

	  if ( sourceReactionNo == lossReactionNo )
	    {
	      sourceReactionNoMultiplier = listOfReactionNoSourcesMultiplier[is];
	      lossReactionNoMultiplier = listOfReactionNoLossesMultiplier[il];

	      if ( sourceReactionNoMultiplier == lossReactionNoMultiplier )
		{ //set both multipliers to zero for later removal
		  listOfReactionNoSourcesMultiplier[is] = 0;
		  listOfReactionNoLossesMultiplier[il] = 0;
		}
	      
	      if ( sourceReactionNoMultiplier > lossReactionNoMultiplier )
		{ //more on source than loss, the subtract multipliers source-loss and set loss multiplier to zero for later removal
		  listOfReactionNoSourcesMultiplier[is] = sourceReactionNoMultiplier - lossReactionNoMultiplier;
		  listOfReactionNoLossesMultiplier[il] = 0;
		}
	      
	      if ( sourceReactionNoMultiplier < lossReactionNoMultiplier )
		{ //more on loss than source, the subtract multipliers loss-source and set source multiplier to zero for later removal
		  listOfReactionNoSourcesMultiplier[is] = 0;
		  listOfReactionNoLossesMultiplier[il] = lossReactionNoMultiplier - sourceReactionNoMultiplier;
		}	      
	    }
	}//il
    }//is


  int listReactionNo[maxNoReactions+1] = {0};
  int listMultiplier[maxNoReactions+1] = {0}; 
  int index = 0;

  /*Tidy up Sources*/
  index = 0;
  for (is=1; is<=NoReactionsSources; is++)
    {
      sourceReactionNoMultiplier = listOfReactionNoSourcesMultiplier[is];
      if (sourceReactionNoMultiplier > 0)
	{
	  index++;
	  listReactionNo[index] = listOfReactionNoSources[is];
	  listMultiplier[index] = listOfReactionNoSourcesMultiplier[is];
	}
    }
  listReactionNo[0] = index;
  listMultiplier[0] = 0;
  
  for (is=0; is<=maxNoReactions; is++)
    {
      listOfReactionNoSources[is] = listReactionNo[is];
      listOfReactionNoSourcesMultiplier[is] = listMultiplier[is];

      listReactionNo[is] = 0;
      listMultiplier[is] = 0;
    }

  /*Tidy up Losses*/
  index = 0;
  for (il=1; il<=NoReactionsLosses; il++)
    {
      lossReactionNoMultiplier = listOfReactionNoLossesMultiplier[il];
      if (lossReactionNoMultiplier > 0)
	{
	  index++;
	  listReactionNo[index] = listOfReactionNoLosses[il];
	  listMultiplier[index] = listOfReactionNoLossesMultiplier[il];
	}
    }
  listReactionNo[0] = index;
  listMultiplier[0] = 0;
  
  for (il=0; il<=maxNoReactions; il++)
    {
      listOfReactionNoLosses[il] = listReactionNo[il];
      listOfReactionNoLossesMultiplier[il] = listMultiplier[il];      
    }

}
