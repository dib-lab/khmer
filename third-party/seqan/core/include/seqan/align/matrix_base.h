// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Simple matrices;  Used in many alignment algorithms.
// ==========================================================================

#ifndef SEQAN_HEADER_MATRIX_BASE_H
#define SEQAN_HEADER_MATRIX_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct NDimensional;


template <typename TValue, unsigned DIMENSION = 0/*typename TSpec = NDimensional*/>
class Matrix;

//////////////////////////////////////////////////////////////////////////////
template <typename T> struct SizeArr_;

template <typename TValue, unsigned DIMENSION>
struct SizeArr_<Matrix<TValue, DIMENSION> >
{
	typedef Matrix<TValue, DIMENSION> TMatrix_;
	typedef typename Size<TMatrix_>::Type TSize_;
	typedef String<TSize_> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
struct Host<Matrix<TValue, DIMENSION> >
{
	typedef String<TValue> Type;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Add more comprehensive documentation!

/*!
 * @class Matrix
 * @headerfile <seqan/align.h>
 * @brief A simple n-dimensional matrix type.
 *
 * @signature template <typename TValue, unsigned DIMENSION = 0>
 *            class Matrix;
 *
 * @tparam TValue    Type of matrix entries.
 * @tparam DIMENSION Dimension of the matrix.  Use 0 for n-dimensional, values &gt; 0 for a matrix with
 *                   <tt>DIMENSION</tt> dimensions.
 */

/**
.Class.Matrix:
..cat:Miscellaneous
..summary:A simple n-dimensional matrix type.
..signature:Matrix<TValue, unsigned DIMENSION = 0>
..param.TValue:Type of matrix entries.
..param.unsigned DIMENSION:The specializing type (0: NDimensional matrix; 2: two dimensional matrix).
..remarks: The following operators and functions are supported: A*B, A*a, A+B,A-B,<<, transpose
..include:seqan/align.h
*/


template <typename TValue>
class Matrix<TValue, 0>
{
//____________________________________________________________________________

public:
	typedef typename Size<Matrix>::Type TSize;
	typedef String<TSize> TSizeArr;
	typedef String<TValue> THost;

	TSizeArr data_lengths;		//Length of every dimension
	TSizeArr data_factors;		//used for positions of dimensions in host ("size of jumps" to get to next entry of specified dimension)

	Holder<THost> data_host;
//____________________________________________________________________________

public:
	Matrix()
	{
		create(data_host);
	}
	Matrix(Matrix const & other_):
		data_lengths(other_.data_lengths),
		data_factors(other_.data_factors),
		data_host(other_.data_host)
	{
	}
	inline Matrix const &
	operator = (Matrix const & other_)
	{
		data_lengths = other_.data_lengths;
		data_factors = other_.data_factors;
		data_host = other_.data_host;

		return *this;
	}
	~Matrix()
	{
	}
//____________________________________________________________________________


//____________________________________________________________________________

	inline TValue &
	operator () (TSize x1, TSize x2)
	{
		return value(*this, x1, x2);
	}
	inline TValue &
	operator () (TSize x1, TSize x2, TSize x3)
	{
		return value(*this, x1, x2, x3);
	}
	inline TValue &
	operator () (TSize x1, TSize x2, TSize x3, TSize x4)
	{
		return value(*this, x1, x2, x3, x4);
	}

//____________________________________________________________________________
};


template <typename TValue>
class Matrix<TValue, 2>
{
//____________________________________________________________________________

public:
	typedef typename Size<Matrix>::Type TSize;
	typedef String<TSize> TSizeArr;
	typedef String<TValue> THost;

	TSizeArr data_lengths;
	TSizeArr data_factors;

	Holder<THost> data_host;


//____________________________________________________________________________

public:
	Matrix()
	{
		create(data_host);

		//setDimension to 2
		resize(data_lengths, 2, 0);
		resize(data_factors, 2, 0);
		data_factors[0] = 1;
	}
	Matrix(Matrix const & other_):
		data_lengths(other_.data_lengths),
		data_factors(other_.data_factors),
		data_host(other_.data_host)
	{
	}
	inline Matrix const &
	operator = (Matrix const & other_)
	{
		data_lengths = other_.data_lengths;
		data_factors = other_.data_factors;
		data_host = other_.data_host;

		return *this;
	}

	~Matrix()
	{
	}
//____________________________________________________________________________


//____________________________________________________________________________

	inline TValue &
	operator () (TSize x1, TSize x2)
	{
		return value(*this, x1, x2);
	}

//____________________________________________________________________________
};

template <typename TValue>
class Matrix<TValue, 3>
{
//____________________________________________________________________________

public:
	typedef typename Size<Matrix>::Type TSize;
	typedef String<TSize> TSizeArr;
	typedef String<TValue> THost;

	TSizeArr data_lengths;
	TSizeArr data_factors;

	Holder<THost> data_host;


//____________________________________________________________________________

public:
	Matrix()
	{
		create(data_host);

		//setDimension to 3
		resize(data_lengths, 3, 0);
		resize(data_factors, 3);
		data_factors[0] = 1;
	}
	Matrix(Matrix const & other_):
		data_lengths(other_.data_lengths),
		data_factors(other_.data_factors),
		data_host(other_.data_host)
	{
	}
	inline Matrix const &
	operator = (Matrix const & other_)
	{
		data_lengths = other_.data_lengths;
		data_factors = other_.data_factors;
		data_host = other_.data_host;

		return *this;
	}

	~Matrix()
	{
	}
//____________________________________________________________________________


//____________________________________________________________________________

	inline TValue &
	operator () (TSize x1, TSize x2, TSize x3)
	{
		return value(*this, x1, x2, x3);
	}

//____________________________________________________________________________
};

template <typename TValue, unsigned DIMENSION>
inline typename SizeArr_<Matrix<TValue, DIMENSION> >::Type &
_dataLengths(Matrix<TValue, DIMENSION> & me)
{
	return me.data_lengths;
}

template <typename TValue, unsigned DIMENSION>
inline typename SizeArr_<Matrix<TValue, DIMENSION> >::Type const &
_dataLengths(Matrix<TValue, DIMENSION> const & me)
{
	return me.data_lengths;
}

template <typename TValue, unsigned DIMENSION>
inline typename SizeArr_<Matrix<TValue, DIMENSION> >::Type &
_dataFactors(Matrix<TValue, DIMENSION> & me)
{
	return me.data_factors;
}

template <typename TValue, unsigned DIMENSION>
inline typename SizeArr_<Matrix<TValue, DIMENSION> >::Type const &
_dataFactors(Matrix<TValue, DIMENSION> const & me)
{
	return me.data_factors;
}

//____________________________________________________________________________


template <typename TValue, unsigned DIMENSION>
inline bool
dependent(Matrix<TValue, DIMENSION> & me)
{
	return dependent(me.data_host);
}

//____________________________________________________________________________

template <typename TValue, unsigned DIMENSION, typename THost>
inline void
setHost(Matrix<TValue, DIMENSION> & me, THost & host_)
{
	setValue(me.data_host, host_);
}

//____________________________________________________________________________


template <typename TValue, unsigned DIMENSION>
inline typename Host<Matrix<TValue, DIMENSION> >::Type &
host(Matrix<TValue, DIMENSION> & me)
{
	return value(me.data_host);
}

template <typename TValue, unsigned DIMENSION>
inline typename Host<Matrix<TValue, DIMENSION> >::Type const &
host(Matrix<TValue, DIMENSION> const & me)
{
	return value(me.data_host);
}

//____________________________________________________________________________

template <typename TValue, unsigned DIMENSION, typename THost>
inline void
assignHost(Matrix<TValue, DIMENSION> & me, THost const & value_)
{
	assignValue(me.data_host, value_);
}

//____________________________________________________________________________

template <typename TValue, unsigned DIMENSION, typename THost>
inline void
moveHost(Matrix<TValue, DIMENSION> & me, THost const & value_)
{
	moveValue(me.data_host, value_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
struct Value< Matrix<TValue, DIMENSION> >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TIteratorSpec>
struct Iterator< Matrix<TValue, DIMENSION>, TIteratorSpec >
{
	typedef Iter<Matrix<TValue, DIMENSION>, PositionIterator> Type;
};

template <typename TValue, unsigned DIMENSION, typename TIteratorSpec>
struct Iterator< Matrix<TValue, DIMENSION> const, TIteratorSpec >
{
	typedef Iter<Matrix<TValue, DIMENSION> const, PositionIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
inline unsigned int
dimension(Matrix<TValue, DIMENSION> const & me)
{
	return length(_dataLengths(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
inline void
setDimension(Matrix<TValue, DIMENSION> & me,
			 unsigned int dim_)
{

	SEQAN_ASSERT_GT(dim_, 0u);
//std::cout<<"\npress enter1\n";
//std::cin.get();
	resize(_dataLengths(me), dim_, 0);

	resize(_dataFactors(me), dim_);
	_dataFactors(me)[0] = 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
inline typename Size<Matrix<TValue, DIMENSION> >::Type
length(Matrix<TValue, DIMENSION> const & me,
	   unsigned int dim_)
{
	return me.data_lengths[dim_];
}

template <typename TValue, unsigned DIMENSION>
inline typename Size<Matrix <TValue, DIMENSION> >::Type
length(Matrix<TValue, DIMENSION> const & me)
{
	return length(host(me));
}

template <typename TValue, unsigned DIMENSION>
inline bool empty(Matrix<TValue, DIMENSION> const & me)
{
	return empty(host(me));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TLength>
inline void
setLength(Matrix<TValue, DIMENSION> & me,
		  unsigned int dim_,
		  TLength length_)
{
	SEQAN_ASSERT_GT(length_, static_cast<TLength>(0));
	SEQAN_ASSERT_LT(dim_, dimension(me));

    typedef typename SizeArr_<Matrix<TValue, DIMENSION> >::TSize_ TSize_;

	_dataLengths(me)[dim_] = static_cast<TSize_>(length_);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
inline void
resize(Matrix<TValue, DIMENSION> & me)
{
	typedef Matrix<TValue, DIMENSION> TMatrix;
	typedef typename Size<TMatrix>::Type TSize;

	unsigned int dimension_ = dimension(me);

	SEQAN_ASSERT_GT(dimension_, 0u);

	TSize factor_ = _dataFactors(me)[0] * length(me, 0);
	for (unsigned int i = 1; (factor_ > 0) && (i < dimension_); ++i)
	{
		_dataFactors(me)[i] = factor_;
		factor_ *= length(me, i);
	}

	if (factor_ > 0)
	{
		resize(host(me), factor_);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TFillValue>
inline void
resize(Matrix<TValue, DIMENSION> & me, TFillValue myValue)	//resize the matrix and fill with value
{
	typedef Matrix<TValue, DIMENSION> TMatrix;
	typedef typename Size<TMatrix>::Type TSize;

	unsigned int dimension_ = dimension(me);

	SEQAN_ASSERT_GT(dimension_, 0u);

	TSize factor_ = _dataFactors(me)[0] * length(me, 0);
	for (unsigned int i = 1; (factor_ > 0) && (i < dimension_); ++i)
	{
		_dataFactors(me)[i] = factor_;
		factor_ *= length(me, i);
	}

	if (factor_ > 0)
		resize(host(me), factor_, myValue);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TPosition>
inline typename Position<Matrix <TValue, DIMENSION> >::Type
nextPosition(Matrix<TValue, DIMENSION> & me,
			 TPosition position_,
			 unsigned int dimension_)
{
	return position_ + _dataFactors(me)[dimension_];
}

template <typename TValue, unsigned DIMENSION, typename TPosition>
inline typename Position<Matrix <TValue, DIMENSION> >::Type
nextPosition(Matrix<TValue, DIMENSION> const & me,
			 TPosition position_,
			 unsigned int dimension_)
{
	return position_ + _dataFactors(me)[dimension_];
}

template <typename TValue, unsigned DIMENSION, typename TPosition>
inline typename Position<Matrix <TValue, DIMENSION> >::Type
previousPosition(Matrix<TValue, DIMENSION> & me,
				 TPosition position_,
				 unsigned int dimension_)
{
	return position_ - _dataFactors(me)[dimension_];
}

template <typename TValue, unsigned DIMENSION, typename TPosition>
inline typename Position<Matrix <TValue, DIMENSION> >::Type
previousPosition(Matrix<TValue, DIMENSION> const & me,
				 TPosition position_,
				 unsigned int dimension_)
{
	return position_ - _dataFactors(me)[dimension_];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TPosition>
inline typename Size< Matrix <TValue, DIMENSION> >::Type
coordinate(Matrix<TValue, DIMENSION> const & me,
		   TPosition position_,
		   unsigned int dimension_)
{
	SEQAN_ASSERT_LT(dimension_, dimension(me));

	if (dimension_ < dimension(me) - 1)
	{
		return (position_ / _dataFactors(me)[dimension_]) % _dataFactors(me)[dimension_ + 1];
	}
	else
	{
		return position_ / _dataFactors(me)[dimension_];
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TTag>
inline typename Iterator<Matrix <TValue, DIMENSION>, Tag<TTag> const>::Type
begin(Matrix<TValue, DIMENSION> & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, DIMENSION>, Tag<TTag> const >::Type(me, 0);
}
template <typename TValue, unsigned DIMENSION, typename TTag>
inline typename Iterator<Matrix <TValue, DIMENSION> const, Tag<TTag> const>::Type
begin(Matrix<TValue, DIMENSION> const & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, DIMENSION> const, Tag<TTag> const >::Type(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TTag>
inline typename Iterator<Matrix <TValue, DIMENSION>, Tag<TTag> const >::Type
end(Matrix<TValue, DIMENSION> & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, DIMENSION>, Tag<TTag> const >::Type(me, length(host(me)));
}
template <typename TValue, unsigned DIMENSION, typename TTag>
inline typename Iterator<Matrix <TValue, DIMENSION> const, Tag<TTag> const >::Type
end(Matrix<TValue, DIMENSION> const & me,
	  Tag<TTag> const)
{
	return typename Iterator<Matrix <TValue, DIMENSION>, Tag<TTag> const >::Type(me, length(host(me)));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TPosition>
inline typename Reference<Matrix<TValue, DIMENSION> >::Type
value(Matrix<TValue, DIMENSION> & me,
	  TPosition position_)
{
	return value(host(me), position_);
}

template <typename TValue, unsigned DIMENSION, typename TPosition>
inline typename Reference<Matrix<TValue, DIMENSION> const>::Type
value(Matrix<TValue, DIMENSION> const & me,
      TPosition position_)
{
    return value(host(me), position_);
}

//____________________________________________________________________________

//two dimensional value access
template <typename TValue, unsigned DIMENSION, typename TOrdinate1, typename TOrdinate2>
inline typename Reference<Matrix<TValue, DIMENSION> >::Type
value(Matrix<TValue, DIMENSION> & me,
	  TOrdinate1 i1,
	  TOrdinate2 i2)
{
	return value(host(me), i1 + i2 * _dataFactors(me)[1]);
}

template <typename TValue, unsigned DIMENSION, typename TOrdinate1, typename TOrdinate2>
inline typename Reference<Matrix<TValue, DIMENSION> const>::Type
value(Matrix<TValue, DIMENSION> const & me,
      TOrdinate1 i1,
      TOrdinate2 i2)
{
    return value(host(me), i1 + i2 * _dataFactors(me)[1]);
}

//____________________________________________________________________________

//3 dimensional value access

template <typename TValue, unsigned DIMENSION, typename TOrdinate1, typename TOrdinate2, typename TOrdinate3>
inline typename Reference<Matrix<TValue, DIMENSION> >::Type
value(Matrix<TValue, DIMENSION> & me,
	  TOrdinate1 i1,
	  TOrdinate2 i2,
	  TOrdinate3 i3)
{
	return value(host(me), i1 + i2 * _dataFactors(me)[1] + i3 * _dataFactors(me)[2]);
}

//____________________________________________________________________________

//4 dimensional value access

template <typename TValue, unsigned DIMENSION, typename TOrdinate1, typename TOrdinate2, typename TOrdinate3, typename TOrdinate4>
inline typename Reference<Matrix<TValue, DIMENSION> >::Type
value(Matrix<TValue, DIMENSION> & me,
	  TOrdinate1 i1,
	  TOrdinate2 i2,
	  TOrdinate3 i3,
	  TOrdinate4 i4)
{
	return value(host(me), i1 + i2 * _dataFactors(me)[1] + i3 * _dataFactors(me)[2] + i4 * _dataFactors(me)[3]);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Iterator: goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
inline void
goNext(Iter<Matrix<TValue, DIMENSION>, PositionIterator> & me,
	   unsigned int dimension_)
{
	setPosition(me, nextPosition(container(me), position(me), dimension_));
}

template <typename TValue, unsigned DIMENSION>
inline void
goNext(Iter<Matrix<TValue, DIMENSION> const, PositionIterator> & me,
	   unsigned int dimension_)
{
	setPosition(me, nextPosition(container(me), position(me), dimension_));
}

template <typename TValue, unsigned DIMENSION>
inline void
goNext(Iter<Matrix<TValue, DIMENSION>, PositionIterator> & me)
{
	goNext(me, 0);
}

template <typename TValue, unsigned DIMENSION>
inline void
goNext(Iter<Matrix<TValue, DIMENSION> const, PositionIterator> & me)
{
	goNext(me, 0);
}

//////////////////////////////////////////////////////////////////////////////
// Iterator: goPrevious
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION>
inline void
goPrevious(Iter< Matrix<TValue, DIMENSION>, PositionIterator > & me,
		   unsigned int dimension_)
{
	setPosition(me, previousPosition(container(me), position(me), dimension_));
}

template <typename TValue, unsigned DIMENSION>
inline void
goPrevious(Iter< Matrix<TValue, DIMENSION> const, PositionIterator > & me,
		   unsigned int dimension_)
{
	setPosition(me, previousPosition(container(me), position(me), dimension_));
}

template <typename TValue, unsigned DIMENSION>
inline void
goPrevious(Iter< Matrix<TValue, DIMENSION>, PositionIterator > & me)
{
	goPrevious(me, 0);
}

template <typename TValue, unsigned DIMENSION>
inline void
goPrevious(Iter< Matrix<TValue, DIMENSION> const, PositionIterator > & me)
{
	goPrevious(me, 0);
}

//////////////////////////////////////////////////////////////////////////////
// goTo
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned DIMENSION, typename TPosition0, typename TPosition1>
inline void
goTo(Iter<Matrix<TValue, DIMENSION>, PositionIterator> & me, TPosition0 pos0, TPosition1 pos1)
{
    setPosition(me, pos0 + pos1 * _dataFactors(container(me))[1]);
}


template <typename TValue, unsigned DIMENSION, typename TPosition0, typename TPosition1>
inline void
goTo(Iter<Matrix<TValue, DIMENSION> const, PositionIterator> & me, TPosition0 pos0, TPosition1 pos1)
{
    setPosition(me, pos0 + pos1 * _dataFactors(container(me))[1]);
}


template <typename TValue, unsigned DIMENSION, typename TPosition0, typename TPosition1, typename TPosition2>
inline void
goTo(Iter<Matrix<TValue, DIMENSION>, PositionIterator> & me, TPosition0 pos0, TPosition1 pos1, TPosition2 pos2)
{
    setPosition(me, pos0 + pos1 * _dataFactors(container(me))[1] + pos2 * _dataFactors(container(me))[2]);
}


template <typename TValue, unsigned DIMENSION, typename TPosition0, typename TPosition1, typename TPosition2>
inline void
goTo(Iter<Matrix<TValue, DIMENSION> const, PositionIterator> & me, TPosition0 pos0, TPosition1 pos1, TPosition2 pos2)
{
    setPosition(me, pos0 + pos1 * _dataFactors(container(me))[1] + pos2 * _dataFactors(container(me))[2]);
}

//////////////////////////////////////////////////////////////////////////////
// Iterator: coordinate

template <typename TValue, unsigned DIMENSION>
inline typename Size< Matrix<TValue, DIMENSION> >::Type
coordinate(Iter<Matrix<TValue, DIMENSION>, PositionIterator > & me,
		   unsigned int dimension_)
{
	return coordinate(container(me), position(me), dimension_);
}

template <typename TValue, unsigned DIMENSION>
inline typename Size< Matrix<TValue, DIMENSION> >::Type
coordinate(Iter<Matrix<TValue, DIMENSION> const, PositionIterator > & me,
		   unsigned int dimension_)
{
	return coordinate(container(me), position(me), dimension_);
}

/*
operator +
Computes the matricial sum between two matrices
..signature:Matrix +(matrix1,matrix2)
..param.matrix1:The first matrix.
...type:Class.Matrix
..param.matrix2:The second matrix.
...type:Class.Matrix
..returns:The sum of the two matrices (another nxm matrix).
..remarks:The number of rows and columns of matrix1 must be equal to the number of rows and columns of matrix2 (length of dimensions for NDimensional matrices)
*/

template <typename TValue,unsigned DIMENSION>
Matrix<TValue,DIMENSION>
operator + (Matrix<TValue,DIMENSION> const & matrix1,Matrix<TValue,DIMENSION> const & matrix2)
{
	//the two matrices must have same dimension
	SEQAN_ASSERT(_dataLengths(matrix1) == _dataLengths(matrix2));

	Matrix<TValue,DIMENSION> result;
	//copy the first matrix
	setDimension(result,length(_dataLengths(matrix1)));
	_dataLengths(result) = _dataLengths(matrix1);
	resize(result);

	//add the matrices
	for(unsigned int i = 0;i< length(host(result));++i)
	{
		value(host(result), i)=value(host(matrix1), i)+value(host(matrix2), i);
	}
	//Return matrix sum
	return result;
}

template <typename TValue,unsigned DIMENSION>
Matrix<TValue,DIMENSION>
operator - (Matrix<TValue,DIMENSION> const & matrix1,Matrix<TValue,DIMENSION> const & matrix2)
{
	//the two matrices must have same dimension
	SEQAN_ASSERT(_dataLengths(matrix1) == _dataLengths(matrix2));

	Matrix<TValue,DIMENSION> result;
	//resize the matrix
	setDimension(result,length(_dataLengths(matrix1)));
	_dataLengths(result) = _dataLengths(matrix1);
	resize(result);

	//subtract the matrices
	for(unsigned int i = 0;i< length(host(result));++i)
	{
		value(host(result), i)=value(host(matrix1), i)-value(host(matrix2), i);
	}
	//Return matrix difference
	return result;
}

template <typename TValue>
Matrix<TValue, 2>
operator * (Matrix<TValue, 2> const & matrix1, Matrix<TValue, 2> const & matrix2)
{
	SEQAN_ASSERT_EQ(length(matrix1,1), length(matrix2,0));

	unsigned int nrow1=length(matrix1,0);
	unsigned int ncol2=length(matrix2,1);
	Matrix<TValue, 2> result;
	//resize the matrix
	setLength(result, 0, nrow1);
	setLength(result, 1, ncol2);
	resize(result,(TValue) 0);

	//Matrix product
	for(unsigned int row = 0; row < nrow1; row++)
	{
		for(unsigned int col = 0; col < ncol2; col++)
		{
			for(unsigned int colRes = 0; colRes < length(matrix1,1); colRes++)
			{
				value(result,row,col)+=	value(host(matrix1), row + colRes * matrix1.data_factors[1])*value(host(matrix2), colRes + col * matrix2.data_factors[1]);
			}
		}
	}
	//return the matrix product
	return result;
}


template <typename TValue>
Matrix<TValue, 2>
operator * (TValue const & scalar, Matrix<TValue, 2> const & matrix)
{
	Matrix<TValue, 2> result;
	result= matrix;
	//scalar multiplication
	for(unsigned int i = 0;i< length(host(result));++i)
	{
		value(host(result), i)*=scalar;
	}
	//return the matrix product
	return result;
}

template <typename TValue>
Matrix<TValue, 2>
operator * (Matrix<TValue, 2> const & matrix, TValue const & scalar)
{
	Matrix<TValue, 2> result;
	result= matrix;
	//scalar multiplication
	for(unsigned int i = 0;i< length(host(result));++i)
	{
		value(host(result), i)*=scalar;
	}
	//return the matrix product
	return result;
}


template <typename TValue, unsigned DIMENSION1, unsigned DIMENSION2>
bool
operator == (Matrix<TValue, DIMENSION1> const & matrix1, Matrix<TValue, DIMENSION2> const & matrix2)
{
	bool result;
	result= (matrix1.data_lengths==matrix2.data_lengths)&&(matrix1.data_factors==matrix2.data_factors)&&(value(matrix1.data_host)==value(matrix2.data_host))&&(DIMENSION1==DIMENSION2);
	return result;
}
/*
.Function.matricialSum:
..summary:Computes the matricial sum between two nxm matrixes
..signature:matricialSum(matrix1,matrix2)
..param.matrix1:The first matrix.
...type:Matrix<TValue, 2>&
..param.matrix2:The second matrix.
...type:Matrix<TValue, 2>&
..returns:The sum of the two matrices (another nxm matrix).
..remarks:The number of rows and columns of matrix1 must be equal to the number of rows and columns of matrix2.
..include:seqan/align.h
*/
/*
template <typename TValue>
Matrix<TValue,2>
matricialSum(Matrix<TValue,2> &matrix1,Matrix<TValue,2> &matrix2)
{
	//the two matrices must have same dimension
	if(length(matrix1,0) != length(matrix2,0)||length(matrix1,1) != length(matrix2,1))
	{
		fprintf(stderr,"Error: The two matrices have different dimensions");
	}


	unsigned int nrow=length(matrix1,0);
	unsigned int ncol=length(matrix1,1);

	Matrix<TValue,2> result;
	//resize the matrix
	setLength(result, 0, nrow);
	setLength(result, 1, ncol);
	resize(result);

	//add the matrices
	for(unsigned int i = 0;i< nrow*ncol;++i)
	{
		value(host(result), i)=value(host(matrix1), i)+value(host(matrix2), i);
	}
	//Return matrix difference
	return result;

}
*/
//////////////////////////////////////////////////////////////////////////////
// _matricialDifference
//////////////////////////////////////////////////////////////////////////////

/*
.Function.matricialDifference:
..summary:Computes the matricial difference between two matrixes
..signature:matricialDifference(matrix1,matrix2)
..param.matrix1:The first matrix.
...type:Matrix<TValue, 2>&
..param.matrix2:The second matrix.
...type:Matrix<TValue, 2>&
..returns:The difference of the two matrices (another matrix).
..remarks:The number of rows and columns of matrix1 must be equal to the number of rows and columns of matrix2.
..include:seqan/align.h
*/
/*
template <typename TValue>
inline Matrix<TValue,2>
matricialDifference(Matrix<TValue,2> & matrix1, Matrix<TValue,2> & matrix2)
{
	//the two matrices must have same dimension
	if(length(matrix1,0) != length(matrix2,0)||length(matrix1,1) != length(matrix2,1))
	{
		fprintf(stderr,"Error: The two matrices have different dimensions");
	}

	unsigned int nrow=length(matrix1,0);
	unsigned int ncol=length(matrix1,1);

	Matrix<TValue,2> result;
	//resize the matrix
	//setDimension(result, 2);
	setLength(result, 0, nrow);
	setLength(result, 1, ncol);
	resize(result);

	//Substract the matrices
	for(unsigned int i1 = 0;i1< nrow;++i1)
		{
			for(unsigned int i2 = 0;i2<ncol;++i2)
			{
				value(host(result), i1 + i2 * _dataFactors(result)[1])=value(host(matrix1), i1 + i2 * _dataFactors(matrix1)[1])-value(host(matrix2), i1 + i2 * _dataFactors(matrix2)[1]);
			}

		}
	//Return matrix difference
	return result;
}
*/
/*
.Function.matricialProduct:
..summary:Computes the matricial product between two 2-dimensional matrixes
..signature:matricialProduct(matrix1,matrix2)
..param.matrix1:The first matrix (mxn).
...type:Matrix<TValue,2>&
..param.matrix2:The second matrix (nxp).
...type:Matrix<TValue,2>&
..returns:The products of the two matrices (another matrix, mxp).
..remarks:The number of columns of matrix1 (left matrix) must be equal to the number of rows of matrix2(right matrix).
..include:seqan/align.h
*/
/*
template <typename TValue>
inline Matrix<TValue, 2>
matricialProduct(Matrix<TValue, 2> &matrix1,
		Matrix<TValue, 2> &matrix2)
{
	//SEQAN_ASSERT_LT(dimension_, dimension(me));
	if(length(matrix1,1) != length(matrix2,0))
	{
		fprintf(stderr,"Error: Number of columns of matrix1 is unequal to number of rows of matrix2");
	}

	unsigned int nrow1=length(matrix1,0);
	unsigned int ncol2=length(matrix2,1);
	Matrix<TValue, 2> result;
	//resize the matrix
	setLength(result, 0, nrow1);
	setLength(result, 1, ncol2);
	resize(result,(TValue) 0);

	//Matrix product
	for(unsigned int row = 0; row < nrow1; row++)
	{
		for(unsigned int col = 0; col < ncol2; col++)
		{
			for(unsigned int colRes = 0; colRes < length(matrix1,1); colRes++)
			{
				value(result,row,col)+=value(matrix1, row,colRes)*value(matrix2,colRes,col);
			}
		}
	}
	//return the matrix product
	return result;
}
*/
// TODO(holtgrew): Should work as the graph-transpose.
/**
.Function.Matrix#transpose
..summary:Transposes matrix
..class:Class.Matrix
..signature:Matrix transpose(matrix)
..param.matrix:The matrix (mxn) to transpose.
...type:Class.Matrix
...remarks: must be of type Matrix<TValue,2> (two dimensional)
..returns:Transposed matrix
..remarks:Only works on two dimensional matrices
..include:seqan/align.h
*/
template <typename TValue>
Matrix<TValue,2>
transpose(Matrix<TValue,2> const & matrix)
{

	unsigned int nrow=length(matrix,0);
	unsigned int ncol=length(matrix,1);

	Matrix<TValue,2> result;
	//resize the matrix
	setLength(result, 0, ncol);
	setLength(result, 1, nrow);
	resize(result);


	for(unsigned int i1 = 0;i1< nrow;++i1)
	{
		for(unsigned int i2 = 0;i2<ncol;++i2)
		{
			value(host(result), i2 + i1 * _dataFactors(result)[1])=value(host(matrix), i1 + i2 * matrix.data_factors[1]);
		}

	}

	//Return transposed matrix
	return result;

}


template < typename TValue >
::std::ostream& operator<<(::std::ostream &out, const Matrix<TValue,2> &matrix)
{
	for(unsigned int i1 = 0;i1< matrix.data_lengths[0];++i1)
	{
			for(unsigned int i2 = 0;i2<(matrix.data_lengths[1]-1);++i2)
			{
				out<<value(host(matrix), i1 + i2 * matrix.data_factors[1])<<"\t";
			}
			//Last line is followd by \n instead of \t
			out<<value(host(matrix), i1 + (matrix.data_lengths[1]-1) *matrix.data_factors[1])<<"\n";
		}

    return out;
}
//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///// READ
/*
 * TODO(goeke) only square matrices of fixed size can be read in...
 */
///////////////////////////////////////////////////////////////
// template < typename TValue >
// void read(FILE *file, Matrix<TValue,2> & matrix)
// {
// 	//unsigned int column_size=3;
// 	unsigned int column_size=pow(4,5);
// 	//read the transition matrix
// 	setLength(matrix, 0, column_size);
// 	setLength(matrix, 1, column_size);
// resize(matrix,0.0);
// 	for(unsigned int row=0; row<column_size; row++)
// 	{
// 		for(unsigned int col=0; col<column_size; col++)
// 		{
// 		  fscanf(file,"%lf ", & value(matrix, row,col));
// 		}
// 		fscanf(file,"\n");
// 	}
// }

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
