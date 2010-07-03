#ifdef LADA_PROPAGATE_STDVECTOR
#  undef LADA_PROPAGATE_STDVECTOR
#endif

// Propagates std::vector member functions, since stl should not be derived from (?).
#define LADA_PROPAGATE_STDVECTOR( propagate_type, propagate_name ) \
  typedef propagate_type :: const_iterator const_iterator; \
  typedef propagate_type :: iterator iterator; \
  typedef propagate_type :: reverse_iterator reverse_iterator; \
  typedef propagate_type :: const_reverse_iterator const_reverse_iterator; \
  typedef propagate_type :: value_type value_type; \
  typedef propagate_type :: reference reference; \
  typedef propagate_type :: const_reference const_reference; \
  typedef propagate_type :: size_type size_type; \
  typedef propagate_type :: difference_type difference_type; \
  const_iterator begin() const { return propagate_name.begin(); } \
  const_iterator end() const { return propagate_name.end(); } \
  iterator begin() { return propagate_name.begin(); } \
  iterator end() { return propagate_name.end(); } \
  const_reverse_iterator rbegin() const { return propagate_name.rbegin(); } \
  const_reverse_iterator rend() const { return propagate_name.rend(); } \
  reverse_iterator rbegin() { return propagate_name.rbegin(); } \
  reverse_iterator rend() { return propagate_name.rend(); } \
  size_type size() const { return propagate_name.size(); } \
  size_type max_size() const { return propagate_name.max_size(); } \
  size_type capacity() const { return propagate_name.capacity(); } \
  bool empty() const { return propagate_name.empty(); } \
  const_reference operator[](size_type _n) const { return propagate_name[_n]; } \
  reference operator[](size_type _n) { return propagate_name[_n]; } \
  reference front() { return propagate_name.front(); } \
  const_reference front() const { return propagate_name.front(); } \
  reference back() { return propagate_name.back(); } \
  const_reference back() const { return propagate_name.back(); } \
  void push_back(value_type const& _t) { propagate_name.push_back(_t); } \
  void pop_back() { propagate_name.pop_back(); } \
  iterator insert(iterator _pos, value_type const& _t) { return propagate_name.insert(_pos, _t); } \
  template<class T_IT> \
    void insert(iterator _pos, T_IT _i1, T_IT _i2) \
      { propagate_name.insert(_pos, _i1, _i2); } \
  iterator erase(iterator _pos) { return propagate_name.erase(_pos); } \
  iterator erase(iterator _i1, iterator _i2) { return propagate_name.erase(_i1, _i2); } \
  void clear() { propagate_name.clear(); } \
  void resize(size_type _n, value_type _t = value_type()) { propagate_name.resize(_n, _t); } \
  void reserve(size_type _n) { propagate_name.reserve(_n); } 
