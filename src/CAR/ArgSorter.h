#ifndef _ArgSorter_
#define _ArgSorter_

// Compare two objects when given pointers to them
template<class Type>
bool compare(const Type* t1, const Type* t2)
{
	return(*(t1) < *(t2));
}

template<class Type>
class ArgSorter
{
	private:
		vector<const Type*> pointers;

	public:
		ArgSorter(const vector<Type>& the_objects)
		:pointers(the_objects.size())
		{
			for(size_t i=0; i<pointers.size(); i++)
				pointers[i] = &(the_objects[i]);
		}

		void sort()
		{
			std::sort(pointers.begin(), pointers.end(), compare<Type>);
		}

		const Type& operator [] (int index) const
		{
			return *(pointers[index]);
		}

		size_t size() const
		{
			return pointers.size();
		}
};

#endif

