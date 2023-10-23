#ifndef BZ_ARRAY_FLUJOS_STENCIL_ET_H
#define BZ_ARRAY_FLUJOS_STENCIL_ET_H

BZ_NAMESPACE(blitz)

template<typename T_ArrayNumtype, int N_rank, typename T_result>
class StencilExpr2 
{
public:
    typedef T_result T_numtype;
    typedef Array<T_ArrayNumtype,N_rank> T_array;
    typedef const T_array& T_ctorArg1;
    typedef int T_ctorArg2;

    enum { numArrayOperands = 2, numIndexPlaceholders = 0,
        rank = N_rank };

    StencilExpr2(const T_array& array1,const T_array& array2)
        : iter1_(array1), iter2_(array2)
    { }

    ~StencilExpr2()
    { }

    // operator* must be declared by subclass
  
    int ascending(int rank)
    { 
		return bounds::compute_ascending(rank, 
		  								 iter1_.ascending(rank), 
									     iter2_.ascending(rank) ); 
	}
 
    int ordering(int rank)
    {
		return bounds::compute_ordering(rank, 
										iter1_.ordering(rank), 
										iter2_.ordering(rank));
	}
 
    int lbound(int rank)
    { 
		return bounds::compute_lbound(rank, 
									  iter1_.lbound(rank), 
									  iter2_.lbound(rank));
	}

    int ubound(int rank)
    {
		return bounds::compute_ubound(rank, 
									  iter1_.ubound(rank), 
									  iter2_.ubound(rank));
	}

    void push(int position)
    { 
		iter1_.push(position);
		iter2_.push(position); 
	}

    void pop(int position)
    {
		iter1_.pop(position); 
		iter2_.pop(position); 
	}

    void advance()
    { 
		iter1_.advance();
		iter2_.advance(); 
	}

    void advance(int n)
    { 
		iter1_.advance(n); 
		iter2_.advance(n); 
	}

    void loadStride(int rank)
    { 
		iter1_.loadStride(rank);
		iter2_.loadStride(rank); 
	}

    bool isUnitStride(int rank) const
    { 
		return iter1_.isUnitStride(rank) 
            && iter2_.isUnitStride(rank);
	}

    void advanceUnitStride()
    { 
		iter1_.advanceUnitStride(); 
		iter2_.advanceUnitStride(); 
	}

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExpr<>::canCollapse()");
		return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank);
    }

    // T_numtype operator[](int i)   -- don't know how to do that.

    // T_numtype fastRead(int i)     -- ditto

    int suggestStride(int rank) const
    { 
		int stride1 = iter1_.suggestStride(rank);
        int stride2 = iter2_.suggestStride(rank);
        return stride1>stride2?stride1:stride2;
	}

    bool isStride(int rank, int stride) const
    { 
		return iter1_.isStride(rank,stride) 
            && iter2_.isStride(rank,stride);
	}

    void prettyPrint(string& str) const
    {
        str += "(stencil)";    // lame, needs work
    }

    void prettyPrint(string& str, prettyPrintFormat& format) const
    {   str += "(stencil)"; }

    template<class T_shape>
    bool shapeCheck(const T_shape& shape)
    { 
		int t1 = iter1_.shapeCheck(shape);
        int t2 = iter2_.shapeCheck(shape);

        return t1 && t2;
	}

    void moveTo(const TinyVector<int,N_rank>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
    }

protected:
    FastArrayIterator<T_ArrayNumtype,N_rank> iter1_, iter2_;
};

template<typename T_ArrayNumtype, int N_rank, typename T_result>
class StencilExpr3 
{
public:
    typedef T_result T_numtype;
    typedef Array<T_ArrayNumtype,N_rank> T_array;
    typedef const T_array& T_ctorArg1;
    typedef int T_ctorArg2;

    enum { numArrayOperands = 3, numIndexPlaceholders = 0,
        rank = N_rank };

    StencilExpr3(const T_array& array1,const T_array& array2,const T_array& array3)
        : iter1_(array1), iter2_(array2), iter3_(array3)
    { }

    ~StencilExpr3()
    { }

    // operator* must be declared by subclass
  
    int ascending(int rank)
    { return bounds::compute_ascending(rank, bounds::compute_ascending(
		rank, iter1_.ascending(rank), iter2_.ascending(rank)),
		iter3_.ascending(rank) ); 
	}
 
    int ordering(int rank)
    {
		return bounds::compute_ordering(rank, bounds::compute_ordering(
          rank, iter1_.ordering(rank), iter2_.ordering(rank)),
		  iter3_.ordering(rank));
	}
 
    int lbound(int rank)
    { 
		return bounds::compute_lbound(rank, bounds::compute_lbound(
          rank, iter1_.lbound(rank), iter2_.lbound(rank)),
		  iter3_.lbound(rank));
	}

    int ubound(int rank)
    {
		return bounds::compute_ubound(rank, bounds::compute_ubound(
          rank, iter1_.ubound(rank), iter2_.ubound(rank)),
		  iter3_.ubound(rank));
	}

    void push(int position)
    { 
		iter1_.push(position);
		iter2_.push(position); 
		iter3_.push(position); 
	}

    void pop(int position)
    {
		iter1_.pop(position); 
		iter2_.pop(position); 
		iter3_.pop(position); 
	}

    void advance()
    { 
		iter1_.advance();
		iter2_.advance(); 
		iter3_.advance(); 
	}

    void advance(int n)
    { 
		iter1_.advance(n); 
		iter2_.advance(n); 
		iter3_.advance(n); 
	}

    void loadStride(int rank)
    { 
		iter1_.loadStride(rank);
		iter2_.loadStride(rank); 
		iter3_.loadStride(rank); 
	}

    bool isUnitStride(int rank) const
    { 
		return iter1_.isUnitStride(rank) 
            && iter2_.isUnitStride(rank)
            && iter3_.isUnitStride(rank);
	}

    void advanceUnitStride()
    { 
		iter1_.advanceUnitStride(); 
		iter2_.advanceUnitStride(); 
		iter3_.advanceUnitStride(); 
	}

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExpr<>::canCollapse()");
		return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank)
			&& iter3_.canCollapse(outerLoopRank, innerLoopRank);
    }

    // T_numtype operator[](int i)   -- don't know how to do that.

    // T_numtype fastRead(int i)     -- ditto

    int suggestStride(int rank) const
    { 
		int stride1 = iter1_.suggestStride(rank);
        int stride2 = iter2_.suggestStride(rank);
        int stride3 = iter3_.suggestStride(rank);
        return stride1>(stride2=(stride2>stride3?stride2:stride3))?stride1:stride2;
	}

    bool isStride(int rank, int stride) const
    { 
		return iter1_.isStride(rank,stride) 
            && iter2_.isStride(rank,stride)
			&& iter3_.isStride(rank,stride);
	}

    void prettyPrint(string& str) const
    {
        str += "(stencil)";    // lame, needs work
    }

    void prettyPrint(string& str, prettyPrintFormat& format) const
    {   str += "(stencil)"; }

    template<class T_shape>
    bool shapeCheck(const T_shape& shape)
    { 
		int t1 = iter1_.shapeCheck(shape);
        int t2 = iter2_.shapeCheck(shape);
        int t3 = iter3_.shapeCheck(shape);

        return t1 && t2 && t3;
	}

    void moveTo(const TinyVector<int,N_rank>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
        iter3_.moveTo(i);
    }

protected:
    FastArrayIterator<T_ArrayNumtype,N_rank> iter1_, iter2_, iter3_;
};

template<typename T_ArrayNumtype, int N_rank, typename T_result>
class StencilExpr4 
{
public:
    typedef T_result T_numtype;
    typedef Array<T_ArrayNumtype,N_rank> T_array;
    typedef const T_array& T_ctorArg1;
    typedef int T_ctorArg2;

    enum { numArrayOperands = 4, numIndexPlaceholders = 0,
        rank = N_rank };

    StencilExpr4(const T_array& array1,const T_array& array2,const T_array& array3,const T_array& array4)
        : iter1_(array1), iter2_(array2), iter3_(array3), iter4_(array4)
    { }

    ~StencilExpr4()
    { }

    // operator* must be declared by subclass
  
    int ascending(int rank)
    { 
		return bounds::compute_ascending(rank,
			bounds::compute_ascending(rank, iter1_.ascending(rank), iter2_.ascending(rank)),
			bounds::compute_ascending(rank, iter3_.ascending(rank), iter4_.ascending(rank))); 
	}
 
    int ordering(int rank)
    {
		return bounds::compute_ordering(rank, 
			bounds::compute_ordering(rank, iter1_.ordering(rank), iter2_.ordering(rank)),
		    bounds::compute_ordering(rank, iter3_.ordering(rank), iter4_.ordering(rank)));
	}
 
    int lbound(int rank)
    { 
		return bounds::compute_lbound(rank, 
			bounds::compute_lbound(rank, iter1_.lbound(rank), iter2_.lbound(rank)),
		    bounds::compute_lbound(rank, iter3_.lbound(rank), iter4_.lbound(rank)));
	}

    int ubound(int rank)
    {
		return bounds::compute_ubound(rank, 
			bounds::compute_ubound(rank, iter1_.ubound(rank), iter2_.ubound(rank)),
		    bounds::compute_ubound(rank, iter3_.ubound(rank), iter4_.ubound(rank)));
	}

    void push(int position)
    { 
		iter1_.push(position);
		iter2_.push(position); 
		iter3_.push(position);
		iter4_.push(position);
	}

    void pop(int position)
    {
		iter1_.pop(position); 
		iter2_.pop(position); 
		iter3_.pop(position); 
		iter4_.pop(position); 
	}

    void advance()
    { 
		iter1_.advance();
		iter2_.advance(); 
		iter3_.advance(); 
		iter4_.advance(); 
	}

    void advance(int n)
    { 
		iter1_.advance(n); 
		iter2_.advance(n); 
		iter3_.advance(n); 
		iter4_.advance(n); 
	}

    void loadStride(int rank)
    { 
		iter1_.loadStride(rank);
		iter2_.loadStride(rank); 
		iter3_.loadStride(rank); 
		iter4_.loadStride(rank); 
	}

    bool isUnitStride(int rank) const
    { 
		return iter1_.isUnitStride(rank) 
            && iter2_.isUnitStride(rank)
            && iter3_.isUnitStride(rank)
            && iter4_.isUnitStride(rank);
	}

    void advanceUnitStride()
    { 
		iter1_.advanceUnitStride(); 
		iter2_.advanceUnitStride(); 
		iter3_.advanceUnitStride(); 
		iter4_.advanceUnitStride(); 
	}

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExpr<>::canCollapse()");
		return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank)
            && iter3_.canCollapse(outerLoopRank, innerLoopRank)
			&& iter4_.canCollapse(outerLoopRank, innerLoopRank);
    }

    // T_numtype operator[](int i)   -- don't know how to do that.

    // T_numtype fastRead(int i)     -- ditto

    int suggestStride(int rank) const
    { 
		int stride1 = iter1_.suggestStride(rank);
        int stride2 = iter2_.suggestStride(rank);
        int stride3 = iter3_.suggestStride(rank);
        int stride4 = iter4_.suggestStride(rank);
		int res = stride1>stride2?stride1:stride2;
		res = stride3>res?stride3:res;
		return stride4>res?stride4:res;
	}

    bool isStride(int rank, int stride) const
    { 
		return iter1_.isStride(rank,stride) 
            && iter2_.isStride(rank,stride)
			&& iter3_.isStride(rank,stride)
			&& iter4_.isStride(rank,stride);
	}

    void prettyPrint(string& str) const
    {
        str += "(stencil)";    // lame, needs work
    }

    void prettyPrint(string& str, prettyPrintFormat& format) const
    {   str += "(stencil)"; }

    template<class T_shape>
    bool shapeCheck(const T_shape& shape)
    { 
		int t1 = iter1_.shapeCheck(shape);
        int t2 = iter2_.shapeCheck(shape);
        int t3 = iter3_.shapeCheck(shape);
        int t4 = iter4_.shapeCheck(shape);

        return t1 && t2 && t3 && t4;
	}

    void moveTo(const TinyVector<int,N_rank>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
        iter3_.moveTo(i);
        iter4_.moveTo(i);
    }

protected:
    FastArrayIterator<T_ArrayNumtype,N_rank> iter1_, iter2_, iter3_, iter4_;
};


template<typename T_ArrayNumtype1, typename T_ArrayNumtype2, int N_rank>
class StencilExpr2Types
{
public:
    typedef T_ArrayNumtype1 T_numtype;
    typedef Array<T_ArrayNumtype1,N_rank> T_array1;
    typedef Array<T_ArrayNumtype2,N_rank> T_array2;
    typedef const T_array1& T_ctorArg1;
    typedef const T_array2& T_ctorArg2;

    enum { numArrayOperands = 2, numIndexPlaceholders = 0,
        rank = N_rank };

    StencilExpr2Types(const T_array1& array1,const T_array2& array2)
        : iter1_(array1), iter2_(array2)
    { }

    ~StencilExpr2Types()
    { }

    // operator* must be declared by subclass
  
    int ascending(int rank)
    { 
		return bounds::compute_ascending(rank, 
		  								 iter1_.ascending(rank), 
									     iter2_.ascending(rank) ); 
	}
 
    int ordering(int rank)
    {
		return bounds::compute_ordering(rank, 
										iter1_.ordering(rank), 
										iter2_.ordering(rank));
	}
 
    int lbound(int rank)
    { 
		return bounds::compute_lbound(rank, 
									  iter1_.lbound(rank), 
									  iter2_.lbound(rank));
	}

    int ubound(int rank)
    {
		return bounds::compute_ubound(rank, 
									  iter1_.ubound(rank), 
									  iter2_.ubound(rank));
	}

    void push(int position)
    { 
		iter1_.push(position);
		iter2_.push(position); 
	}

    void pop(int position)
    {
		iter1_.pop(position); 
		iter2_.pop(position); 
	}

    void advance()
    { 
		iter1_.advance();
		iter2_.advance(); 
	}

    void advance(int n)
    { 
		iter1_.advance(n); 
		iter2_.advance(n); 
	}

    void loadStride(int rank)
    { 
		iter1_.loadStride(rank);
		iter2_.loadStride(rank); 
	}

    bool isUnitStride(int rank) const
    { 
		return iter1_.isUnitStride(rank) 
            && iter2_.isUnitStride(rank);
	}

    void advanceUnitStride()
    { 
		iter1_.advanceUnitStride(); 
		iter2_.advanceUnitStride(); 
	}

    bool canCollapse(int outerLoopRank, int innerLoopRank) const
    {
        // BZ_DEBUG_MESSAGE("_bz_ArrayExpr<>::canCollapse()");
		return iter1_.canCollapse(outerLoopRank, innerLoopRank)
            && iter2_.canCollapse(outerLoopRank, innerLoopRank);
    }

    // T_numtype operator[](int i)   -- don't know how to do that.

    // T_numtype fastRead(int i)     -- ditto

    int suggestStride(int rank) const
    { 
		int stride1 = iter1_.suggestStride(rank);
        int stride2 = iter2_.suggestStride(rank);
        return stride1>stride2?stride1:stride2;
	}

    bool isStride(int rank, int stride) const
    { 
		return iter1_.isStride(rank,stride) 
            && iter2_.isStride(rank,stride);
	}

    void prettyPrint(string& str) const
    {
        str += "(stencil)";    // lame, needs work
    }

    void prettyPrint(string& str, prettyPrintFormat& format) const
    {   str += "(stencil)"; }

    template<class T_shape>
    bool shapeCheck(const T_shape& shape)
    { 
		int t1 = iter1_.shapeCheck(shape);
        int t2 = iter2_.shapeCheck(shape);

        return t1 && t2;
	}

    void moveTo(const TinyVector<int,N_rank>& i)
    {
        iter1_.moveTo(i);
        iter2_.moveTo(i);
    }

protected:
    FastArrayIterator<T_ArrayNumtype1,N_rank> iter1_;
    FastArrayIterator<T_ArrayNumtype2,N_rank> iter2_;
};


#define BZ_ET_STENCIL2(name,result) \
template<class P_numtype, int N_rank> \
class name ## _et : public StencilExpr2<P_numtype,N_rank,result>, \
  public ETBase<name ## _et<P_numtype,N_rank> > \
 { \
private: \
    typedef StencilExpr2<P_numtype,N_rank,result> T_base; \
    using T_base::iter1_; \
    using T_base::iter2_; \
public: \
    name ## _et(const Array<P_numtype,N_rank>& A,const Array<P_numtype,N_rank>& B) \
        : StencilExpr2<P_numtype,N_rank,result>(A,B) \
    { } \
    result operator*() \
	{ return name(iter1_,iter2_); } \
    result operator()(const TinyVector<int,N_rank>& a) \
    { iter1_.moveTo(a); iter2_.moveTo(a); return name(iter1_,iter2_); } \
    result fastRead(int i) \
    { \
      const P_numtype* tmp1 = iter1_.data(); \
      const P_numtype* tmp2 = iter2_.data(); \
      iter1_._bz_setData(tmp1 + i); \
      iter2_._bz_setData(tmp2 + i); \
      P_numtype r = name(iter1_,iter2_); \
      iter1_._bz_setData(tmp1); \
      iter2_._bz_setData(tmp2); \
      return r; \
    } \
}; \
template<class P_numtype, int N_rank> \
inline _bz_ArrayExpr<name ## _et<P_numtype, N_rank> > \
name(Array<P_numtype,N_rank>& A,Array<P_numtype,N_rank>& B) \
{ \
    return _bz_ArrayExpr<name ## _et<P_numtype, N_rank> >(A,B); \
}

#define BZ_ET_STENCIL3(name,result) \
template<class P_numtype, int N_rank> \
class name ## _et : public StencilExpr3<P_numtype,N_rank,result>, \
  public ETBase<name ## _et<P_numtype,N_rank> > \
 { \
private: \
    typedef StencilExpr3<P_numtype,N_rank,result> T_base; \
    using T_base::iter1_; \
    using T_base::iter2_; \
    using T_base::iter3_; \
public: \
    name ## _et(const Array<P_numtype,N_rank>& A,const Array<P_numtype,N_rank>& B,const Array<P_numtype,N_rank>& C) \
        : StencilExpr3<P_numtype,N_rank,result>(A,B,C) \
    { } \
    result operator*() \
	{ return name(iter1_,iter2_,iter3_); } \
    result operator()(const TinyVector<int,N_rank>& a) \
    { iter1_.moveTo(a); iter2_.moveTo(a); iter3_.moveTo(a); return name(iter1_,iter2_,iter3_); } \
    result fastRead(int i) \
    { \
      const P_numtype* tmp1 = iter1_.data(); \
      const P_numtype* tmp2 = iter2_.data(); \
      const P_numtype* tmp3 = iter3_.data(); \
      iter1_._bz_setData(tmp1 + i); \
      iter2_._bz_setData(tmp2 + i); \
      iter3_._bz_setData(tmp3 + i); \
      P_numtype r = name(iter1_,iter2_,iter3_); \
      iter1_._bz_setData(tmp1); \
      iter2_._bz_setData(tmp2); \
      iter3_._bz_setData(tmp3); \
      return r; \
    } \
}; \
template<class P_numtype, int N_rank> \
inline _bz_ArrayExpr<name ## _et<P_numtype, N_rank> > \
name(Array<P_numtype,N_rank>& A,Array<P_numtype,N_rank>& B,Array<P_numtype,N_rank>& C) \
{ \
    return _bz_ArrayExpr<name ## _et<P_numtype, N_rank> >(A,B,C); \
}

#define BZ_ET_STENCIL4(name,result) \
template<class P_numtype, int N_rank> \
class name ## _et : public StencilExpr4<P_numtype,N_rank,result>, \
  public ETBase<name ## _et<P_numtype,N_rank> > \
 { \
private: \
    typedef StencilExpr4<P_numtype,N_rank,result> T_base; \
    using T_base::iter1_; \
    using T_base::iter2_; \
    using T_base::iter3_; \
    using T_base::iter4_; \
public: \
    name ## _et(const Array<P_numtype,N_rank>& A,const Array<P_numtype,N_rank>& B,const Array<P_numtype,N_rank>& C,const Array<P_numtype,N_rank>& D) \
        : StencilExpr4<P_numtype,N_rank,result>(A,B,C,D) \
    { } \
    result operator*() \
	{ return name(iter1_,iter2_,iter3_,iter4_); } \
    result operator()(const TinyVector<int,N_rank>& a) \
    { iter1_.moveTo(a); iter2_.moveTo(a); iter3_.moveTo(a); iter4_.moveTo(a); return name(iter1_,iter2_,iter3_,iter4_); } \
    result fastRead(int i) \
    { \
      const P_numtype* tmp1 = iter1_.data(); \
      const P_numtype* tmp2 = iter2_.data(); \
      const P_numtype* tmp3 = iter3_.data(); \
      const P_numtype* tmp4 = iter4_.data(); \
      iter1_._bz_setData(tmp1 + i); \
      iter2_._bz_setData(tmp2 + i); \
      iter3_._bz_setData(tmp3 + i); \
      iter4_._bz_setData(tmp4 + i); \
      P_numtype r = name(iter1_,iter2_,iter3_,iter4_); \
      iter1_._bz_setData(tmp1); \
      iter2_._bz_setData(tmp2); \
      iter3_._bz_setData(tmp3); \
      iter4_._bz_setData(tmp4); \
      return r; \
    } \
}; \
template<class P_numtype, int N_rank> \
inline _bz_ArrayExpr<name ## _et<P_numtype, N_rank> > \
name(Array<P_numtype,N_rank>& A,Array<P_numtype,N_rank>& B,Array<P_numtype,N_rank>& C,Array<P_numtype,N_rank>& D) \
{ \
    return _bz_ArrayExpr<name ## _et<P_numtype, N_rank> >(A,B,C,D); \
}

#define BZ_ET_STENCIL_DIFF2(name) \
template<class P_numtype, int N_rank> \
class name ## _et : public StencilExpr<P_numtype,N_rank,P_numtype>, \
  public ETBase<name ## _et<P_numtype,N_rank> > \
 { \
private: \
    typedef StencilExpr<P_numtype,N_rank,P_numtype> T_base; \
    using T_base::iter_; \
public: \
    name ## _et(const Array<P_numtype,N_rank>& A, int dim1, int dim2) \
        : StencilExpr<P_numtype,N_rank,P_numtype>(A), \
          dim1_(dim1), dim2_(dim2) \
    { } \
    P_numtype operator*() \
    { return name(iter_,dim1_,dim2_); } \
    P_numtype operator()(const TinyVector<int,N_rank>& a) \
    { iter_.moveTo(a); return name(iter_,dim1_,dim2_); } \
    P_numtype fastRead(int i) \
    { \
      const P_numtype* tmp = iter_.data(); \
      iter_._bz_setData(tmp + i); \
      P_numtype r = name(iter_,dim1_,dim2_); \
      iter_._bz_setData(tmp); \
      return r; \
    } \
private: \
    int dim1_, dim2_; \
}; \
template<class P_numtype, int N_rank> \
inline _bz_ArrayExpr<name ## _et<P_numtype, N_rank> > \
name(Array<P_numtype,N_rank>& A, int dim1, int dim2) \
{ \
    return _bz_ArrayExpr<name ## _et<P_numtype, N_rank> >(A,dim1,dim2); \
}

// return type	: TinyVector
// arguments	: Array< TinyVector< Pixel, rank > , rank >, Array< Pixel, rank >

#define BZ_ET_STENCILV2(name,rank) \
template<class P_numtype1, class P_numtype2, int N_rank> \
class name ## _et : public StencilExpr2Types<P_numtype1,P_numtype2,N_rank>, \
  public ETBase<name ## _et<P_numtype1,P_numtype2,N_rank> > \
 { \
private: \
    typedef StencilExpr2Types<P_numtype1,P_numtype2,N_rank> T_base; \
    using T_base::iter1_; \
    using T_base::iter2_; \
public: \
    typedef P_numtype1 result; \
    name ## _et(const Array<P_numtype1,N_rank>& A,const Array<P_numtype2,N_rank>& B, int nsize) \
        : StencilExpr2Types<P_numtype1,P_numtype2,N_rank>(A,B), \
          nsize_(nsize) \
	{ } \
    result operator*() \
    { return name(iter1_,iter2_,nsize_); } \
    result operator()(const TinyVector<int,N_rank>& a) \
    { iter1_.moveTo(a); iter2_.moveTo(a); return name(iter1_,iter2_,nsize_); } \
    result fastRead(int i) \
    { \
      const P_numtype1* tmp1 = iter1_.data(); \
      const P_numtype2* tmp2 = iter2_.data(); \
      iter1_._bz_setData(tmp1 + i); \
      iter2_._bz_setData(tmp2 + i); \
      result r = name(iter1_,iter2_,nsize_); \
      iter1_._bz_setData(tmp1); \
      iter2_._bz_setData(tmp2); \
      return r; \
    } \
private: \
    int nsize_; \
}; \
template<class P_numtype1, class P_numtype2, int N_rank> \
inline _bz_ArrayExpr<name ## _et<P_numtype1, P_numtype2, N_rank> > \
name(Array<P_numtype1,N_rank>& A, Array<P_numtype2,N_rank>& B,int nsize) \
{ \
    return _bz_ArrayExpr< name ## _et<P_numtype1,P_numtype2,N_rank> >(A,B,nsize); \
}

BZ_ET_STENCIL_DIFF2(mixed22)
BZ_ET_STENCIL_DIFF2(mixed22n)

BZ_NAMESPACE_END

#endif // BZ_ARRAY_FLUJOS_STENCIL_ET_H
