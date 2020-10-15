#pragma once

#include <iterator>
namespace Foam
{

//- Convert an iterator to a reverse iterator
template <class InputIt>
std::reverse_iterator<InputIt> make_reverse(InputIt i)
{
    return std::reverse_iterator<InputIt>(i);
}

//- Compute the number of sequential elements for which the predicate is true
template <class InputIt, class UnaryPredicate>
typename std::iterator_traits<InputIt>::difference_type
count_while(InputIt first, InputIt last, UnaryPredicate predicate)
{

    typename std::iterator_traits<InputIt>::difference_type ret = 0;
    for(; first != last; ++first)
    {
        if(!predicate(*first))
        {
            break;
        }
        ret++;
    }
    return ret;
}

} // namespace Foam