#ifndef INCLUDE_
#define INCLUDE_

#include <cassert>
#include <cstring>
#include <experimental/iterator>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <vector>
#include <math.h>
#include <omp.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fftw3.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace cbclib {

namespace py = pybind11;

template <typename Container, typename Shape, typename = std::enable_if_t<
    std::is_rvalue_reference_v<Container &&> && std::is_integral_v<typename std::remove_cvref_t<Shape>::value_type>
>>
inline py::array_t<typename Container::value_type> as_pyarray(Container && seq, Shape && shape)
{
    Container * seq_ptr = new Container(std::move(seq));
    auto capsule = py::capsule(seq_ptr, [](void * p) {delete reinterpret_cast<Container *>(p);});
    return py::array(std::forward<Shape>(shape),  // shape of array
                     seq_ptr->data(),  // c-style contiguous strides for Container
                     capsule           // numpy array references this parent
    );
}

template <typename Container, typename = std::enable_if_t<std::is_rvalue_reference_v<Container &&>>>
inline py::array_t<typename Container::value_type> as_pyarray(Container && seq)
{
    Container * seq_ptr = new Container(std::move(seq));
    auto capsule = py::capsule(seq_ptr, [](void * p) {delete reinterpret_cast<Container *>(p);});
    return py::array(seq_ptr->size(),  // shape of array
                     seq_ptr->data(),  // c-style contiguous strides for Container
                     capsule           // numpy array references this parent
    );
}

template <typename Container, typename Shape, typename = std::enable_if_t<
    std::is_integral_v<typename std::remove_cvref_t<Shape>::value_type>
>>
inline py::array_t<typename Container::value_type> to_pyarray(const Container & seq, Shape && shape)
{
    return py::array(std::forward<Shape>(shape), seq.data());
}

template <typename Container>
inline py::array_t<typename Container::value_type> to_pyarray(const Container & seq)
{
    return py::array(seq.size(), seq.data());
}

template <typename T, typename = void>
struct is_input_iterator : std::false_type {};

template <typename T>
struct is_input_iterator<T,
    std::void_t<decltype(*std::declval<T &>()), decltype(++std::declval<T &>())>
> : std::true_type {};

template <typename T>
constexpr bool is_input_iterator_v = is_input_iterator<T>::value;

namespace detail {

struct Constants
{
    // 1 / sqrt(2 * pi)
    static constexpr double M_1_SQRT2PI = 0.3989422804014327;
};

template <typename T>
class any_container
{
protected:
    std::vector<T> vec;

public:
    any_container() = default;

    template <typename It, typename = std::enable_if_t<is_input_iterator<It>::value>>
    any_container(It first, It last) : vec(first, last) {}

    template <typename Container,
        typename = std::enable_if_t<
            std::is_convertible_v<decltype(*std::begin(std::declval<const Container &>())), T>
        >
    >
    any_container(const Container & c) : any_container(std::begin(c), std::end(c)) {}

    // initializer_list's aren't deducible, so don't get matched by the above template;
    // we need this to explicitly allow implicit conversion from one:
    template <typename TIn, typename = std::enable_if_t<std::is_convertible_v<TIn, T>>>
    any_container(const std::initializer_list<TIn> & c) : any_container(c.begin(), c.end()) {}

    any_container(std::vector<T> && v) : vec(std::move(v)) {}

    // Moves the vector out of an rvalue any_container
    operator std::vector<T> && () && { return std::move(this->vec); }

    std::vector<T> & operator*() {return this->vec;}
    const std::vector<T> & operator*() const {return this->vec;}

    std::vector<T> * operator->() {return &(this->vec);}
    const std::vector<T> * operator->() const {return &(this->vec);}
};

}

template <typename Container>
void check_dimensions(const std::string & name, ssize_t axis, const Container & shape) {}

template <typename Container, typename... Ix>
void check_dimensions(const std::string & name, ssize_t axis, const Container & shape, ssize_t i, Ix... index)
{
    if (axis < 0)
    {
        auto text = name + " has the wrong number of dimensions: " + std::to_string(shape.size()) +
                    " < " + std::to_string(shape.size() - axis);
        throw std::invalid_argument(text);
    }
    if (axis >= static_cast<ssize_t>(shape.size()))
    {
        auto text = name + " has the wrong number of dimensions: " + std::to_string(shape.size()) +
                    " < " + std::to_string(axis + 1);
        throw std::invalid_argument(text);
    }
    if (i != static_cast<ssize_t>(shape[axis]))
    {
        auto text = name + " has an incompatible shape at axis " + std::to_string(i) + ": " +
                    std::to_string(shape[axis]) + " != " + std::to_string(i);
        throw std::invalid_argument(text);
    }
    check_dimensions(name, axis + 1, shape, index...);
}

template <typename Container, typename = std::enable_if_t<std::is_integral_v<typename Container::value_type>>>
void fail_container_check(std::string msg, const Container & shape)
{
    std::ostringstream oss;
    std::copy(shape.begin(), shape.end(), std::experimental::make_ostream_joiner(oss, ", "));
    throw std::invalid_argument(msg + ": {" + oss.str() + "}");
}

template <typename ForwardIt1, typename ForwardIt2>
void check_equal(const std::string & msg, ForwardIt1 first1, ForwardIt1 last1, ForwardIt2 first2, ForwardIt2 last2)
{
    if (!std::equal(first1, last1, first2))
    {
        std::ostringstream oss1, oss2;
        std::copy(first1, last1, std::experimental::make_ostream_joiner(oss1, ", "));
        std::copy(first2, last2, std::experimental::make_ostream_joiner(oss2, ", "));
        throw std::invalid_argument(msg + ": {" + oss1.str() + "}, {" + oss2.str() + "}");
    }
}

template <typename T, int ExtraFlags>
void fill_array(py::array_t<T, ExtraFlags> & arr, T fill_value)
{
    auto fill = py::array_t<T, ExtraFlags>(py::ssize_t(1));
    fill.mutable_at(0) = fill_value;
    PyArray_CopyInto(reinterpret_cast<PyArrayObject *>(arr.ptr()), reinterpret_cast<PyArrayObject *>(fill.ptr()));
}

template <typename ForwardIt, typename V, int ExtraFlags>
void check_optional(const std::string & name, ForwardIt first, ForwardIt last,
                    std::optional<py::array_t<V, ExtraFlags>> & opt, V fill_value)
{
    if (!opt)
    {
        opt = py::array_t<V, ExtraFlags>(std::vector(first, last));
        fill_array(opt.value(), fill_value);
    }
    py::buffer_info obuf = opt.value().request();
    check_equal(name + " and inp arrays must have identical shapes",
                obuf.shape.begin(), obuf.shape.end(), first, last);
}

template <typename T>
class sequence : public detail::any_container<T>
{
public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    using detail::any_container<T>::any_container;

    template <typename U, typename = std::enable_if_t<std::is_convertible_v<U, T>>>
    sequence(U value, size_t length = 1)
    {
        std::fill_n(std::back_inserter(this->vec), length, value);
    }

    template <typename Container, typename = std::enable_if_t<std::is_convertible_v<typename Container::value_type, T>>>
    sequence(const Container & vec, size_t length)
    {
        if (vec.size() < length)
            throw std::invalid_argument("rank of vector (" + std::to_string(vec.size()) +
                                        ") is less than the required length (" + std::to_string(length) + ")");
        std::copy_n(vec.begin(), length, std::back_inserter(this->vec));
    }

    size_t size() const {return this->vec.size();}
    sequence & unwrap(T max)
    {
        for (size_t i = 0; i < this->size(); i++)
        {
            this->vec[i] = (this->vec[i] >= 0) ? this->vec[i] : max + this->vec[i];
            if (this->vec[i] >= max)
                throw std::invalid_argument("axis " + std::to_string(this->vec[i]) +
                                            " is out of bounds (ndim = " + std::to_string(max) + ")");
        }
        return *this;
    }

    template <class Array, typename = std::enable_if_t<std::is_base_of_v<py::array, std::remove_cvref_t<Array>>>>
    Array && swap_axes(Array && arr) const
    {
        size_t counter = 0;
        for (py::ssize_t i = 0; i < arr.ndim(); i++)
        {
            if (std::find(this->vec.begin(), this->vec.end(), i) == this->vec.end())
            {
                auto obj = reinterpret_cast<PyArrayObject *>(arr.release().ptr());
                arr = py::reinterpret_steal<std::remove_cvref_t<Array>>(PyArray_SwapAxes(obj, counter++, i));
            }
        }
        return std::forward<Array>(arr);
    }

    template <class Array, typename V = std::remove_cvref_t<Array>::value_type, typename = std::enable_if_t<
        std::is_same_v<py::array_t<V>, std::remove_cvref_t<Array>>
    >>
    Array && swap_axes_back(Array && arr) const
    {
        size_t counter = arr.ndim() - this->size();
        for (py::ssize_t i = arr.ndim() - 1; i >= 0; i--)
        {
            if (std::find(this->vec.begin(), this->vec.end(), i) == this->vec.end())
            {
                auto obj = reinterpret_cast<PyArrayObject *>(arr.release().ptr());
                arr = py::reinterpret_steal<std::remove_cvref_t<Array>>(PyArray_SwapAxes(obj, --counter, i));
            }
        }
        return std::forward<Array>(arr);
    }

    T & operator[] (size_t index) {return this->vec[index];}
    const T & operator[] (size_t index) const {return this->vec[index];}

    iterator begin() {return this->vec.begin();}
    const_iterator begin() const {return this->vec.begin();}
    iterator end() {return this->vec.end();}
    const_iterator end() const {return this->vec.end();}
};

template <typename F, typename = std::enable_if_t<std::is_floating_point<F>::value>>
bool isclose(F a, F b, F atol = F(1e-8), F rtol = F(1e-5))
{
    if (fabs(a - b) <= atol + rtol * std::fmax(fabs(a), fabs(b))) return true;
    return false;
}

inline void * import_numpy() {import_array(); return NULL;}

class thread_exception
{
    std::exception_ptr ptr;
    std::mutex lock;

public:

    thread_exception() : ptr(nullptr) {}

    void rethrow()
    {
        if (this->ptr) std::rethrow_exception(this->ptr);
    }

    void capture_exception()
    {
        std::unique_lock<std::mutex> guard(this->lock);
        this->ptr = std::current_exception();
    }

    template <typename Function, typename... Args>
    void run(Function && f, Args &&... params)
    {
        try
        {
            std::forward<Function>(f)(std::forward<Args>(params)...);
        }
        catch (...)
        {
            capture_exception();
        }
    }
};

}

#endif
