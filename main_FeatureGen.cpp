#include  <iostream>
#include  "stdio.h"
#include  <cmath>
#include  "getCommonSeq/ComSeq.h"
using namespace std;

int main( int argc , char *argv[] )
{

     ClusterToFeature Test;

     if(argv[1]==NULL)
     {
          printf("FeatureGen ToFeature\n");
          printf("FeatureGen Decipher\n");
          printf("FeatureGen Statistics\n");
          printf("FeatureGen CorStatistics\n");
          printf("FeatureGen LogNormalized\n");
          return 0;
     }



     if( strcmp(argv[1],"ToFeature")==0)
     {


        FILE * ip =  fopen(argv[2],"r");
        FILE * op =  fopen(argv[3],"w+");

        if(ip==NULL||op==NULL||argc<5)
        {
            printf("FeatureGen ToFeature   InputFile OutputFIle AverageCorrelation AverageComponentRatio\n");
            return 0;
        }
        Test.tolerance = 0.01;
        Test.CorThreshold = atof(argv[4]);
        Test.CompRThreshold = atof(argv[5]);
        Test.ClusterRead(ip);
        Test.OutPutFeature(op);
     }

     if( strcmp(argv[1],"Decipher")==0)
     {


        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op =  fopen(argv[4],"w+");

         if(ip==NULL||op==NULL||argc<4)
         {
            printf("FeatureGen Decipher FeatureFile   InputFile OutputFile\n");
            return 0;
         }

        Test.ReadFeature(ip);
        Test.SampleRead(ip1);
        Test.SampleDecipher(op);

     }

     if( strcmp(argv[1],"FeatureToNetwork")==0)
     {

        FILE * ip =  fopen(argv[2],"r");
        FILE * op =  fopen(argv[3],"w");
        FILE * op1 =  fopen(argv[4],"w");
        std::map < std::string ,int > GeneName;

        Test.ReadFeature(ip);
        int NODESIZE =0;
        for(int i=0;i<Test.FeatureM.size();i++)
        {

           NODESIZE = NODESIZE + Test.FeatureM[i].n_rows;
        }

        int NodeIndex =0;
        for(int i=0;i< Test.FeatureM.size();i++)
        {
          int geneID;

           for( int j=0;j<Test.FeatureNameArray[i].size();j++)
           {

            if ( GeneName.find( Test.FeatureNameArray[i][j].c_str()) ==  GeneName.end() )
            {

               GeneName.insert(std::make_pair(Test.FeatureNameArray[i][j].c_str(),GeneName.size()+ 1.0 + NODESIZE));
               geneID = GeneName[Test.FeatureNameArray[i][j].c_str()];

            }else
            {
               geneID = GeneName[Test.FeatureNameArray[i][j].c_str()];
            }

            for(int k=0;k<Test.FeatureM[i].n_rows;k++)
            {
              fprintf(op,"%d    ",geneID);
              fprintf(op,"%d    ",NodeIndex+k+1);
              fprintf(op,"%lf\n",Test.FeatureM[i](k,j));
            }

          }
          NodeIndex = NodeIndex + Test.FeatureM[i].n_rows;
        }


         std::map < std::string,int >::iterator iter;
         fprintf(op1,"NODESIZE:");
         fprintf(op1,"%d\n",NODESIZE);
         for (iter = GeneName.begin(); iter != GeneName.end(); ++iter ) {

              string TempCh = iter->first;
              int ID = iter->second;
              fprintf(op1,"%s   ",TempCh.c_str());
              fprintf(op1,"%d\n",ID);
         }

     }

     if( strcmp(argv[1],"Statistics")==0)
     {
        char buffer[1024];
        std::vector < std::vector < unsigned int > > Index;
        std::vector < std::vector < double > > Value;


        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op = fopen(argv[4],"w+");

         if(ip==NULL||op==NULL||argc<4)
         {
            printf("FeatureGen Statistics InputFile InputFileGroupTag OutputFile\n");
            return 0;

         }

        std::vector < double > CorSum;
        std::vector < double > SquareSum;

        Test.ReadSpareM(Index,Value,ip);

        std::vector < unsigned int > Flag;
        Flag.resize(Index.size());

        for(int i=0;i<Index.size();i++ )
        {
            fscanf(ip1,"\t");
            fscanf(ip1,"%d",&Flag[i]);
            fscanf(ip1,"\n");
        }

        Test.StatisticSpareMGroup(Index,Value,Flag,CorSum,SquareSum);

        for(int i=0;i<CorSum.size();i++)
        {
            fprintf(op,"%lf", CorSum[i]);
            fprintf(op,"\t");
            fprintf(op,"%lf", SquareSum[i]);
            fprintf(op,"\n");
        }
     }

     
   if( strcmp(argv[1],"LogNormalized") ==0 )
   {

       FILE *ip = fopen(argv[2],"r"); 
       FILE *op = fopen(argv[3],"w+");


       if(ip==NULL||op==NULL||argc<3)
       {
	      printf("FeatureGen LogNormalized InputFile OutputFile\n");
           return 0;
       }

       Test.SampleRead(ip); 
       Test.Normalized();
       Test.Log();
       Test.SampleWrite(op);
      
   }
     
     
     
      if( strcmp(argv[1],"CorStatistics")==0)
     {
        char buffer[1024];
        std::vector < std::vector < unsigned int > > Index;
        std::vector < std::vector < double > > Value;


        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op = fopen(argv[4],"w+");


         if(ip==NULL||op==NULL||argc<4)
         {
            printf("FeatureGen CorStatistics InputFile InputFileNumericalTag OutputFile\n");
            return 0;
         }



        std::vector < double > CorSum;

        Test.ReadSpareM(Index,Value,ip);

        std::vector < double > Flag;
        Flag.resize(Index.size());

        for(int i=0;i<Index.size();i++ )
        {
            fscanf(ip1,"\t");
            fscanf(ip1,"%lf",&Flag[i]);
            fscanf(ip1,"\n");
        }

        Test.StatisticSpareCor(Index,Value,Flag,CorSum);

        for(int i=0;i<CorSum.size();i++)
        {
            fprintf(op,"%lf",CorSum[i]);
            fprintf(op,"\n");
        }
     }


    if( strcmp(argv[1],"NetworkFilter")==0)
    {

        FILE * ip =  fopen(argv[2],"r");
        FILE * ip1 = fopen(argv[3],"r");
        FILE * op = fopen(argv[4],"w+");


       Test.ReadNetWork(ip1);
       Test.ReadFeature(ip);
       Test.NetWorkFilteFeature();
       Test.OutPutFeatureS(op);
    }
/*
   if( strcmp(argv[1],"NetworkKernal") ==0)
   {
       FILE * ip =  fopen(argv[2],"r");
       FILE * ip1 = fopen(argv[3],"r");
       FILE * op = fopen(argv[4],"w+");

       Test.SampleRead(ip);
       Test.ReadNetWork(ip1);
       Test.NetWorkKernal();
       Test.OutPutNetWorkKernalM(op);
   }

*/
   if( strcmp(argv[1],"NetworkEnvironment") ==0)
   {
       FILE * ip =  fopen(argv[2],"r");
       FILE * ip1 = fopen(argv[3],"r");
       FILE * ip2 = fopen(argv[4],"r");
       FILE * op = fopen(argv[5],"w+");
       char buffer[1024];
       Test.target="T";
       Test.SampleRead(ip);
       Test.CellType.resize(Test.Sample.n_rows);
       for(int i=0;i<Test.CellType.size();i++ )
       {
            fscanf(ip1,"\t");
            fscanf(ip1,"\n");
            fscanf(ip1,"%s",buffer);
            Test.CellType[i].assign(buffer);
            fscanf(ip1,"\n");
       }
       Test.ReadNetWork(ip2);
       Test.NetEnvironment();
       Test.OutPutNetWorkKernalM(op);

       //to do
   }




   return 0;
}




