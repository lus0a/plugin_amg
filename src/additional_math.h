namespace ug{
template<typename TMatrix, typename TVector>
double ANorm(const TMatrix &A, TVector &v)
{
#ifdef UG_PARALLEL
	if(A.has_storage_type(PST_ADDITIVE))
		v.change_storage_type(PST_CONSISTENT);
	else
		v.change_storage_type(PST_ADDITIVE);
#endif
	THROW_IF_NOT_EQUAL(v.size(), A.num_cols());
	double d = 0;
	typename TVector::value_type b;
	for(size_t i=0; i<A.num_rows(); i++)
	{
		b=0;
		A.mat_mult_add_row(i, b, 1.0, v);
		d += VecProd(b, v[i]);
	}

#ifdef UG_PARALLEL
	return sqrt(v.layouts()->proc_comm().allreduce(d, PCL_RO_SUM));
#else
	return sqrt(d);
#endif
}

template<typename TMatrix, typename TVector>
double DNorm(const TMatrix &A, TVector &v)
{
#ifdef UG_PARALLEL
	if(A.has_storage_type(PST_ADDITIVE))
		v.change_storage_type(PST_CONSISTENT);
	else
		v.change_storage_type(PST_ADDITIVE);
#endif
	THROW_IF_NOT_EQUAL(v.size(), A.num_cols());
	double d = 0;
	for(size_t i=0; i<A.num_rows(); i++)
		d += VecProd(A(i, i)*v[i], v[i]);
#ifdef UG_PARALLEL
	return sqrt(v.layouts()->proc_comm().allreduce(d, PCL_RO_SUM));
#else
	return sqrt(d);
#endif
}


template<typename TMatrix, typename TVector>
double DApply(const TMatrix &A, TVector &v)
{
#ifdef UG_PARALLEL
	if(A.has_storage_type(PST_ADDITIVE))
		v.change_storage_type(PST_CONSISTENT);
	else
		v.change_storage_type(PST_ADDITIVE);
#endif
	THROW_IF_NOT_EQUAL(v.size(), A.num_cols());
	for(size_t i=0; i<A.num_rows(); i++)
		v[i] = A(i, i)*v[i];
}




}
