#include "GraphEmbedding.h"


void  laplaceNormalized(sp_mat & AdjacencyM,sp_mat & LaplaceM )
{

     std:vector <double> D;

     D.resize(AdjacencyM.n_rows);
     LaplaceM.resize( AdjacencyM.n_rows, AdjacencyM.n_cols);
     for(int i=0;i<AdjacencyM.n_rows;i++)
     {
       D[i]=0;
     }

     for ( int i=0;i<AdjacencyM.n_rows;i++)
      {
        sp_mat::row_iterator it = AdjacencyM.begin_row(i);
        sp_mat::row_iterator it_end = AdjacencyM.end_row(i);

        for(;it !=it_end;++it)
	    {
          D[i] = D[i]+(*it);
		}
	  }

      for ( int i=0;i<AdjacencyM.n_rows;i++)
	  {
	     D[i] = sqrt(D[i]);
	  }

	  for ( int i=0;i<AdjacencyM.n_rows;i++)
      {
	    sp_mat::row_iterator   it = AdjacencyM.begin_row(i);
		sp_mat::row_iterator   it_end = AdjacencyM.end_row(i);

		for(;it !=it_end;++it)
		{
			int row;
			int col;

			row= it.row();
			col= it.col();
			double value = (*it);
            LaplaceM.at(row,col)=-value/(D[row]*D[col]);
		}
	  }
	  for ( int i=0;i<AdjacencyM.n_rows;i++)
      {
          LaplaceM.at(i,i) = 1+LaplaceM.at(i,i);
      }
}


void  laplace(sp_mat & AdjacencyM,sp_mat & LaplaceM )
{

    std:vector <double> D;

    D.resize(AdjacencyM.n_rows);
    LaplaceM.resize( AdjacencyM.n_rows, AdjacencyM.n_cols);
    for(int i=0;i<AdjacencyM.n_rows;i++)
    {
        D[i]=0;
    }

    for ( int i=0;i<AdjacencyM.n_rows;i++)
    {

       sp_mat::row_iterator  it = AdjacencyM.begin_row(i);
       sp_mat::row_iterator  it_end = AdjacencyM.end_row(i);

        for(;it !=it_end;++it)
	    {
          D[i] = D[i]+(*it);
		}
    }
     for ( int i=0;i<AdjacencyM.n_rows;i++)
     {
	    sp_mat::row_iterator it = AdjacencyM.begin_row(i);
		sp_mat::row_iterator it_end = AdjacencyM.end_row(i);

		for(;it !=it_end;++it)
		{
			int row;
			int col;

			row= it.row();
			col= it.col();
			double value = (*it);
            LaplaceM.at(row,col) = -value;
		}
	  }
     for ( int i=0;i<AdjacencyM.n_rows;i++ )
     {
         LaplaceM.at(i,i) = D[i]+LaplaceM.at(i,i);
     }
}


void  GraphEncode(sp_mat & LaplaceM)
{
		  vec eigval;
		  mat eigvec;
          double x = 0.5;

		  eigs_sym(eigval,eigvec,LaplaceM,LaplaceM.n_rows);
          double T12 = boost::math::chebyshev_t(12, x);

}












