#ifndef ARRAY2D_HH
#define ARRAY2D_HH

#include <cassert>
#include <algorithm>

/**
 * @brief Container class for a basic 2d array
 * @tparam T Datatype to be stored
 */
template <typename T>
struct Array2d
{
    /// Default constructor for an empty, zero-size array
    Array2d() : m_data(nullptr), m_size1(0), m_size2(0), m_owndata(false) { }

    /// Constructs a new array and allocates its memory
    Array2d(size_t size1, size_t size2) :
        m_size1(size1), m_size2(size2), m_owndata(true)
    {
        m_data = new T[m_size1*m_size2];
    }

    /// Basic destructor to clean up heap-allocated memory
    virtual ~Array2d()
    {
        if (m_owndata) delete[] m_data;
    }

    /// Disabled copy constructor
    Array2d(Array2d const& other) = delete;

    /// Disabled copy assignment
    Array2d& operator=(Array2d const& other) = delete;

    /// Move constructor (no copies)
    Array2d(Array2d&& other) : Array2d()
    {
        std::swap(m_data, other.m_data);
        std::swap(m_size1, other.m_size1);
        std::swap(m_size2, other.m_size2);
        std::swap(m_owndata, other.m_owndata);
    }

    /// Setter for element (i,j) of the array
    T& operator()(size_t i, size_t j)
    {
        assert((i < m_size1) && (j < m_size2));
        return m_data[i*m_size2+j];
    }

    /// Getter for element (i,j) of the array
    T operator()(size_t i, size_t j) const
    {
        assert((i < m_size1) && (j < m_size2));
        return m_data[i*m_size2+j];
    }

    /// Getter for the raw data buffer
    T *data() { return m_data; }

    /// Const getter for the raw data buffer
    const T *data() const { return m_data; }

    protected:
        T *m_data;          ///< Raw data pointer
        size_t m_size1;     ///< Array size on dimension 1
        size_t m_size2;     ///< Array size on dimension 2
        bool m_owndata;     ///< Flag for data ownership
};

#endif      // ARRAY2D_HH

// vim: set ft=cpp.doxygen:
