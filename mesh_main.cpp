#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <omp.h>

using namespace std;

#ifndef SAFE_MODE
#define SAFE_MODE
#endif

#ifdef SAFE_MODE
#define SAFE_ASSERT(X,...) if(!(X)) exit(printf(__VA_ARGS__))
                                    // throw?
#else
#define SAFE_ASSERT(X,...)
#endif

struct Node {
    double x;
    double y;
    Node(double _x, double _y){
        x = _x;
        y = _y;
    }
};

typedef int tInt; // index

struct Edge {
    tInt first, second;
    Edge(tInt _first, tInt _second){
        first = _first;
        second = _second;
    }
    Edge() {
        first = 0;
        second = 0;
    }
};

struct Boundary {
    Edge edge;
    int direction;
    Boundary(tInt _first, tInt _second, int _direction) {
        edge.first = _first;
        edge.second = _second;
        direction = _direction;
    }
};

enum tIdxType {
    IDX_N   = 0,  // nodes
    IDX_E   = 1,  // elements
    IDX_S   = 2  // segments
};

template<tIdxType T>
class tIdx {
private:
    tInt i;
    inline tIdx(const tInt j) : i(j) {}

public:
    inline tIdx() : i(-1) {}

    inline void Set(int j) { i = j; }

    inline tInt idx() const { return i; }

    inline tIdx<T>& operator=(const tInt j) { i = j; return *this; }
    inline tIdx<T>& operator=(const size_t j) { i = (tInt)j; return *this; }
    inline tIdx<T>& operator=(const tIdx<T>& j){ i = j.i; return *this; }

    inline tIdx<T>& operator++() { ++i; return *this; }
    inline tIdx<T>& operator--() { --i; return *this; }

    inline tIdx<T> operator++(int) { return tIdx<T>(i++); }
    inline tIdx<T> operator--(int) { return tIdx<T>(i--); }

    inline tIdx<T> operator+(tIdx<T> j) { return tIdx<T>(i + j.i); }
    inline tIdx<T> operator-(tIdx<T> j) { return tIdx<T>(i - j.i); }

    inline tIdx<T> operator+(tInt j) { return tIdx<T>(i + j); }
    inline tIdx<T> operator-(tInt j) { return tIdx<T>(i - j); }
    inline tIdx<T> operator*(tInt j) { return tIdx<T>(i * j); }
    inline tIdx<T> operator/(tInt j) { return tIdx<T>(i / j); }
    inline tIdx<T> operator%(tInt j) { return tIdx<T>(i % j); }

    inline bool operator<(tInt j) const { return i < j; }
    inline bool operator>(tInt j) const { return i > j; }
    inline bool operator<=(tInt j) const { return i <= j; }
    inline bool operator>=(tInt j) const { return i >= j; }
    inline bool operator==(tInt j) const { return i == j; }
    inline bool operator!=(tInt j) const { return i != j; }

    inline bool operator<(tIdx<T> j) const { return i < j.i; }
    inline bool operator>(tIdx<T> j) const { return i > j.i; }
    inline bool operator<=(tIdx<T> j) const { return i <= j.i; }
    inline bool operator>=(tIdx<T> j) const { return i >= j.i; }
    inline bool operator==(tIdx<T> j) const { return i == j.i; }
    inline bool operator!=(tIdx<T> j) const { return i != j.i; }

    inline bool operator!() { return !i; }

    template<tIdxType TT>
    tIdx(const tIdx<TT> &J); // íå ïîíÿòíî íóæíî ëè
};

template<typename T>
class tBlock {
private:
    T *p;
    const char *name; // pointer to parent object's name
    int N;

    inline tBlock() : N(0), p(nullptr), name(nullptr) {}
    inline tBlock<T>& operator=(const tBlock<T>& B) {
        if (this != &B) { N = B.N; p = B.v; name = B.name; }
        return *this;
    }

public:
    inline tBlock(const tBlock<T> &B) : N(B.N), p(B.p), name(B.name) {}
    inline tBlock(T* ptr, int n, const char* label = nullptr) : p(ptr), N(n), name(label) {}

    inline ~tBlock() {}

    inline void Reset() { N = 0; p = nullptr; name = nullptr; }

    inline operator T* () { return p; }
    inline operator const T* () const { return p; }

    inline int Size() const { return N; }
    inline T* Body() const { return p; }

    inline bool Allocated() const { return (N > 0) && (p != nullptr); }

    inline T& operator[](int i) {
        SAFE_ASSERT((p && i >= 0 && i < N), "tBlock: wrong index [%d] size=%d (%s)\n", i, N, name ? name : "");
        return p[i];
    }
    inline const T& operator[](int i) const {
        SAFE_ASSERT((p && i >= 0 && i < N), "tBlock: wrong index [%d] size=%d (%s)\n", i, N, name ? name : "");
        return p[i];
    }
};

// char* name needed?
template <tIdxType V, typename T> // list for every V type (tIdxType) to it's T - type of kept value (tIdx<tIdxType>)
struct Topo {
    vector<tInt> IA;
    vector<T> JA;
    char * name ;

    void clear() {
        IA.clear();
        JA.clear();
        name = nullptr;
    }

   /* Topo(string name_in) {
        name = name_in;
        IA = vector<tInt>(0);
        JA = vector<T>(0);
    }
*/
    int GetBlockSize(tIdx<V> i) {
        SAFE_ASSERT((i >= 0 && i <= (tInt) IA.size()), "Topo: wrong index [%d] for block in (%s)\n", i.idx(), name);
        if (i == (tInt) IA.size() - 1) {
            return (tInt) JA.size() - IA[i.idx()];
        }
        return IA[i.idx() + 1] - IA[i.idx()];
    }

    inline tBlock<T> operator[](const tIdx<V> id) {
        SAFE_ASSERT((id >= 0 && id < IA.size()), "Topo: wrong index [%d] blocks_count=%d (%s)\n", id.idx(), (tInt) IA.size(), name);
        return tBlock<T>(&JA[IA[id.idx()]], GetBlockSize(id));
    }

    T get(tIdx<V> i, int j) {
        return (*this)[i][j];
    }

    void add(vector<T> &a) {
        IA.push_back((tInt) JA.size());
        for (tInt i = 0; i < (tInt) a.size(); ++i) {
            JA.push_back(a[i]);
        }
    }

    int elem_number() {
        return (int) IA.size();
    }

    size_t GetOccupiedMemory() {

        return IA.capacity() * sizeof(tInt) + JA.capacity() * sizeof(T) + ((name)?(size_t)strlen(name):size_t(0));
    }
};


//ìîæåò îí õðàíèò òáëîê??
template <tIdxType V, typename T> // for each V its T's: T - type of kept value, V - idxType of kept value
struct tBlockVector {
    vector<T> JA;
    int M;//ðàçìåð áëîêà
    string name;

    void clear() {
      JA.clear();
    }

    tBlockVector(int m_input, string nm) {
      M = m_input;
      name = nm;
      JA = vector<T>(0);
    }

    void add(vector<T> &a) {
        SAFE_ASSERT((a.size() == M), "tBlockVector: wrong size of input array (%s)", name.c_str());
        for (int i = 0; i < a.size(); ++i) {
            JA.push_back(a[i]);
        }
    }

    int elem_number() {
        return (JA.size() / M);
    }

    //ïðàâèëüíî ëè
    //âàðèàíòû äåëàåì ñòðóêòóðó ãäå õðàíÿòñÿ ññûëêè íà òáëîêè
    //ïîïðàâèòü ñåéô àññåðò
    inline tBlock<T> operator[](const tIdx<V> id) {
        SAFE_ASSERT((id >= 0 && id < elem_number()), "Topo: wrong index [%d] blocks_count=%d (%s)", id.idx(), elem_number(), name.c_str());
        return tBlock<T>(&JA[M*id.idx()], GetBlockSize(id));
    }

};


template <tIdxType V, typename T>
struct tVector {
    vector<T> VeC;//ÈËÈ vector<tIdx<V>>; - åñëè template <tIdxType V> è äàëüøå ïîïðàâèòü
    string name;


    tVector(int n = 0, string nm = "") {
        name = nm;
        VeC = vector<T>(n);
    }

    //ìåíÿåì ðàçìåð âåêòîðà, ìîæåò realloc, ìîæåò äîñòàòî÷íî êîíñòðóêòîðà
    void set_size(int n){
        VeC = vector<T>(n); // ìîæåò resize()
    }

    void add(T t){
        VeC.push_back(t);
    }

    inline T& operator[] (const tIdx<V> in){
        SAFE_ASSERT((in.idx() >= 0 && in.idx() < VeC.size()), "tVector: wrong index [%d] vector size=%d (%s)", in.idx(), VeC.size(), name.c_str());
        return VeC[in.idx()];
    }

    int elem_number() {
        return VeC.size();
    }

    void clear() {
        VeC.clear();
    }
};


//tIdx <IDX_N> ln, gn
//VeC<DX_N, tIdx<IDX_N>> testvector
// testvector[ln] = gn

template <tIdxType V, typename T>
struct tMap {
    int Nglob, Nloc;
    map<tIdx<V>, T> MaP;
    T minusone;
    string name;


    //êàê ñäåëàòü âûäåëåíèå ïî ôèêñèðîâàííîìó ðàçìåðó?
    //íåâîçìîæíî íå çíàÿ êëþ÷à
    /*åñòü òàêîé ïðèìåð
    std::map<keytype, valuetype,="" std::less<keytype="">, MyAllocator> myMap;

    myMap.get_allocator().reserve( nodeSize * numberOfNodes );
    </keytype,>

    èëè ÷åðåç allocator
    */

    void add(tIdx<V> iglob, T iloc) {
      SAFE_ASSERT(iglob.idx() >= 0 && iglob.idx() < Nglob, "tMap: wrong index [%d] size of global range=%d (%s)", iglob.idx(), Nglob, name.c_str());
      SAFE_ASSERT(iloc >= 0 && iloc < Nloc, "tMap: wrong index [%d] size of local range=%d (%s)", iloc.idx(), Nloc, name.c_str());

      MaP.insert(make_pair(iglob, iloc));
    }

    tMap() : Nglob(0), Nloc(0), name(""), MaP() {}

    tMap(int N1, int N2, string name_in = "") : MaP() {
      Nglob = N1;
      Nloc = N2;
      name = name_in;
    }

    void set_size(int n){
      // íè÷åãî íå äåëàåì?
    }

    int elem_number() {
      return MaP.size();
    }

    void clear() {
      MaP.clear();
      Nglob = Nloc = 0;
    }

    inline T& operator[] (const tIdx<V> iglob) {
      typename map<tIdx<V>, T>::iterator iter;
      iter = MaP.find(iglob);
      if (iter == MaP.end())
        return minusone;
      return iter->second;
    }
};


template <tIdxType V, typename T>
struct tUnorderedMap {
    int Nglob, Nloc;
    unordered_map<tIdx<V>, T> UnordMaP;
    T minusone;
    string name;

    void add(tIdx<V> iglob, T iloc) {
      SAFE_ASSERT(iglob.idx() >= 0 && iglob.idx() < Nglob, "tMap: wrong index [%d] size of global range=%d (%s)", iglob.idx(), Nglob, name.c_str());
      SAFE_ASSERT(iloc >= 0 && iloc < Nloc, "tMap: wrong index [%d] size of local range=%d (%s)", iloc.idx(), Nloc, name.c_str());

      UnordMaP.insert(make_pair(iglob, iloc));
    }

    tUnorderedMap() : Nglob(0), Nloc(0), name(""), UnordMaP() {}

    tUnorderedMap(int N1, int N2, string name_in = "") : UnordMaP() {
      Nglob = N1;
      Nloc = N2;
      name = name_in;
    }

//êàê íàïèñàòü? êàê ñðàçó ñîçäàòü îïðåäåë¸ííîãî ðàçìåðà
    void set_size(int n) {
    }

    int elem_number() {
      return UnordMaP.size();
    }

    void clear () {
      UnordMaP.clear();
      Nglob = Nloc = 0;
    }

    inline T& operator[] (const tIdx<V> iglob) {
      typename map<tIdx<V>, T>::iterator iter;
      /*â ïðèìåðå òóò ñòîÿëî auto*/iter = UnordMaP.find(iglob);
      if (iter == UnordMaP.end())
        return minusone;
      return iter->second;
    }
};

struct Mesh {
public:
    vector<Node> nodes;
    Topo< IDX_E, tIdx<IDX_N> > T; // T == ENadj
    Topo< IDX_N, tIdx<IDX_E> > RT; // RT == NEadj
    Topo< IDX_N, tIdx<IDX_N> > NNadj;
    Topo< IDX_N, tIdx<IDX_S> > NSadj;
    Topo< IDX_S, tIdx<IDX_E> > SEadj;
    Topo< IDX_E, tIdx<IDX_S> > ESadj;
    vector<Boundary> boundary;
    vector<Edge> edges;
    Mesh *child, *father;
    map<tInt, tInt> nodeConnectionWithChild;
    vector<tInt> nodeConnectionWithFather;

    size_t GetOccupiedMemory() {
        return nodes.capacity() * sizeof(Node) + T.GetOccupiedMemory() + RT.GetOccupiedMemory()
            + NNadj.GetOccupiedMemory() + NSadj.GetOccupiedMemory() + SEadj.GetOccupiedMemory()
            + ESadj.GetOccupiedMemory() + boundary.capacity() * sizeof(Boundary)
            + edges.capacity() * sizeof(Edge) + nodeConnectionWithFather.capacity() * sizeof(tInt)
            + sizeof(nodeConnectionWithChild) + nodeConnectionWithChild.size() * sizeof(pair<tInt, tInt>);
    }

    void print_all_debug() {
        cout << "NODES" << endl;
        for (tInt i = 0; i < (tInt) nodes.size(); ++i) {
            cout << nodes[i].x << " " << nodes[i].y << endl;
        }
        cout << "ELEMENTS->NODES" << endl;
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                cout << T[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "NODES->ELEMENTS" << endl;
        tIdx<IDX_N> ne_beg, ne_end;
        ne_beg = 0;
        ne_end = RT.elem_number();
        for (tIdx<IDX_N> i = ne_beg; i < ne_end; ++i) {
            for (int j = 0; j < RT.GetBlockSize(i); ++j) {
                cout << RT[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "Node to Node adjacency" << endl;
        tIdx<IDX_N> nn_beg, nn_end;
        nn_beg = 0;
        nn_end = NNadj.elem_number();
        for (tIdx<IDX_N> i = nn_beg; i < nn_end; ++i) {
            for (int j = 0; j < NNadj.GetBlockSize(i); ++j) {
                cout << NNadj[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "Node to Segment adjacency" << endl;
        tIdx<IDX_N> ns_beg, ns_end;
        ns_beg = 0;
        ns_end = NSadj.elem_number();
        for (tIdx<IDX_N> i = ns_beg; i < ns_end; ++i) {
            for (int j = 0; j < NSadj.GetBlockSize(i); ++j) {
                cout << NSadj[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "EDGES" << endl;
        for (tInt i = 0; i < (tInt) edges.size(); ++i) {
            cout << edges[i].first << " " << edges[i].second << endl;
        }
        cout << "BOUNDARY" << endl;
        for (tInt i = 0; i < (tInt) boundary.size(); ++i) {
            cout << boundary[i].edge.first << " " << boundary[i].edge.second << " " << boundary[i].direction << endl;
        }
    }

/// gens

    void gen_inverse_topo() {
        RT.clear();
        vector<tInt> e(nodes.size());
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                ++e[T[i][j].idx()];
            }
        }
        for (tInt i = 0; i < (tInt) e.size(); ++i) {
            tInt curSize = (tInt) RT.JA.size();
            RT.IA.push_back(curSize);
            tIdx<IDX_E> elem;
            elem = e[i];
            RT.JA.push_back(elem);
            RT.JA.resize(curSize + e[i]);
        }

        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                tIdx<IDX_N> elem = T.get(i, j);
                --RT.JA[RT.IA[elem.idx()]];
                RT.JA[ RT.IA[elem.idx()] + RT.JA[RT.IA[elem.idx()]].idx() ] = i;
            }
        }
    }

    void gen_NNadj() {
        NNadj.clear();
        tIdx<IDX_N> ne_beg, ne_end;
        ne_beg = 0;
        ne_end = RT.elem_number();
        for (tIdx<IDX_N> i = ne_beg; i < ne_end; ++i) {
            unordered_set<tInt> un;
            vector< tIdx<IDX_N> > v;
            for (int j = 0; j < RT.GetBlockSize(i); ++j) {
                tIdx<IDX_E> elem = RT.get(i, j);
                int elem_size = T.GetBlockSize(elem);
                tInt in_elem = -1;
                for (int k = 0; k < elem_size; ++k) {
                    if (i == T.get(elem, k)) {
                        in_elem = k;
                        break;
                    }
                }
                for (int k = 0; k < elem_size; ++k) {
                    tIdx<IDX_N> nn = T.get(elem, k);
                    if ((k != in_elem) && (elem_size == 3 || ((k % 3 == 0) ^ (in_elem % 3 == 0)))) {
                        if (un.find(nn.idx()) == un.end()) {
                            un.insert(nn.idx());
                            v.push_back(nn);
                        }
                    }
                }
            }
            NNadj.add(v);
        }
    }

    void gen_NSadj() {
        NSadj.clear();
        vector<tInt> e(nodes.size(), 0);
        tIdx<IDX_N> ns_beg, ns_end;
        ns_beg = 0;
        ns_end = NSadj.elem_number();
        for (tInt i = 0; i < (tInt) edges.size(); ++i) {
            ++e[edges[i].first];
            ++e[edges[i].second];
        }

        for (tInt i = 0; i < (tInt) e.size(); ++i) {
            tInt curSize = (tInt) NSadj.JA.size();
            NSadj.IA.push_back(curSize);
            tIdx<IDX_S> elem;
            elem = e[i];
            NSadj.JA.push_back(elem);
            NSadj.JA.resize(curSize + e[i]);
        }

        for (tInt i = 0; i < (tInt) edges.size(); ++i) {
            tInt elem = edges[i].first;
            --NSadj.JA[NSadj.IA[elem]];
            NSadj.JA[ NSadj.IA[elem] + NSadj.JA[NSadj.IA[elem]].idx() ] = i;
            elem = edges[i].second;
            --NSadj.JA[NSadj.IA[elem]];
            NSadj.JA[ NSadj.IA[elem] + NSadj.JA[NSadj.IA[elem]].idx() ] = i;
        }
    }

    void gen_edges() {
        edges.clear();
        tIdx<IDX_N> nn_beg, nn_end;
        nn_beg = 0;
        nn_end = NNadj.elem_number();
        for (tIdx<IDX_N> i = nn_beg; i < nn_end; ++i) {
            for (int j = 0; j < NNadj.GetBlockSize(i); ++j) {
                if (i < NNadj.get(i, j)) {
                    edges.push_back(Edge(i.idx(), NNadj.get(i, j).idx()));
                }
            }
        }
    }

    void gen_ESadj() {
        ESadj.clear();
        vector< tIdx<IDX_S> > ans;
        unordered_set<tInt> v;
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            v.clear();
            ans.clear();
            for (int j = 0; j < T.GetBlockSize(i); ++j) { // äîáàâèì âñå óçëû ýëåìåíòà â ñåò
                v.insert(T.get(i, j).idx());
            }
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                tIdx<IDX_N> cur = T.get(i, j); // ïîëó÷èëè íîìåð òåêóùåãî óçëà ýëåìåíòà
                for (int k = 0; k < NSadj.GetBlockSize(cur); ++k) {
                    tIdx<IDX_S> c_edge = NSadj.get(cur, k); // íîìåð òåêóùåãî ðåáðà
                    if (v.count(edges[c_edge.idx()].first) && v.count(edges[c_edge.idx()].second)) {
                        ans.push_back(c_edge);
                    }
                }
                v.erase(cur.idx());
            }
            ESadj.add(ans);
        }
    }

    void gen_SEadj() {
        SEadj.clear();
        SEadj.IA.resize(edges.size());
        for (tInt i = 0; i < (tInt) edges.size(); ++i) {
            SEadj.IA[i] = 2 * i;
        }
        vector< tIdx<IDX_E> > ja(2 * edges.size());
        tIdx<IDX_E> es_beg, es_end;
        es_beg = 0;
        es_end = ESadj.elem_number();
        for (tIdx<IDX_E> i = es_beg; i < es_end; ++i) { // i-ûé ýëåìåíò
            for (int j = 0; j < ESadj.GetBlockSize(i); ++j) { // ïî åãî ðåáðàì èäåì
                ja[SEadj.IA[ESadj.get(i, j).idx()]++] = i;
            }
        }
        for (tInt i = 0; i < (tInt) edges.size(); ++i) {
            tInt n_size = (tInt) SEadj.JA.size();
            for (tInt j = 2 * i; j < SEadj.IA[i]; ++j) {
                SEadj.JA.push_back(ja[j]);
            }
            SEadj.IA[i] = n_size;
        }
    }

/// saves

    void save_coord(int process_number) {
        ofstream out(string("coordinate") + to_string(process_number) + ".msh");
        for (tInt i = 0; i < (tInt) nodes.size(); ++i) {
            out << nodes[i].x << " " << nodes[i].y << endl;
        }
        out.close();
    }

    void save_back_topo(int process_number) {
        ofstream out(string("backtopo") + to_string(process_number) + ".msh");
        tIdx<IDX_N> ne_beg, ne_end;
        ne_beg = 0;
        ne_end = RT.elem_number();
        for (tIdx<IDX_N> i = ne_beg; i < ne_end; ++i) {
            out << RT.GetBlockSize(i);
            for (int j = 0; j < RT.GetBlockSize(i); ++j) {
                out << " " << RT.get(i, j).idx();
            }
            out << endl;
        }
        out.close();
    }

    void save_mesh(int process_number) {
        ofstream out(string("mesh") + to_string(process_number) + ".txt");
        out << "Nn " << nodes.size() << endl << "Nt " << T.elem_number() << endl << "NFaceBC " << boundary.size() << endl << "NumCoords 2" << endl;
        out.close();
    }

    void save_nbf(int process_number) {
        ofstream out(string("bctopo") + to_string(process_number) + ".msh");
        for (tInt i = 0; i < (tInt) boundary.size(); ++i) {
            out << 2 << " " << boundary[i].edge.first << " " << boundary[i].edge.second << " " << boundary[i].direction << endl;
        }
        out.close();
    }

    void save_topo(int process_number) {
        ofstream out(string("topo") + to_string(process_number) + ".msh");
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            out << T.GetBlockSize(i);
            for (int j = 0; j < 3; ++j) {
                out << " " << T.get(i, j).idx();
                if (T.GetBlockSize(i) == 4 && j == 1) {
                    out << " " << T.get(i, 3).idx();
                }
            }
            out << endl;
        }
        out.close();
    }

    void save_NNadj(int process_number) {
        ofstream out(string("NNadj") + to_string(process_number) + ".msh");
        tIdx<IDX_N> nn_beg, nn_end;
        nn_beg = 0;
        nn_end = NNadj.elem_number();
        for (tIdx<IDX_N> i = nn_beg; i < nn_end; ++i) {
            out << NNadj.GetBlockSize(i);
            for (int j = 0; j < NNadj.GetBlockSize(i); ++j) {
                out << " " << NNadj[i][j].idx();
            }
            out << endl;
        }
        out.close();
    }

    void save_ESadj(int process_number) {
        ofstream out(string("ESadj") + to_string(process_number) + ".msh");
        tIdx<IDX_E> es_beg, es_end;
        es_beg = 0;
        es_end = ESadj.elem_number();
        for (tIdx<IDX_E> i = es_beg; i < es_end; ++i) {
            out << ESadj.GetBlockSize(i);
            for (int j = 0; j < ESadj.GetBlockSize(i); ++j) {
                out << " " << ESadj.get(i, j).idx();
            }
            out << endl;
        }
        out.close();
    }

    void save_SEadj(int process_number) {
        ofstream out(string("SEadj") + to_string(process_number) + ".msh");
        tIdx<IDX_S> se_beg, se_end;
        se_beg = 0;
        se_end = SEadj.elem_number();
        for (tIdx<IDX_S> i = se_beg; i < se_end; ++i) { // ó ðåáðà ëèáî 2 ýëåìåíòà ëèáî 1
            out << SEadj.GetBlockSize(i);
            for (int j = 0; j < SEadj.GetBlockSize(i); ++j) {
                out << " " << SEadj.get(i, j).idx();
            }
            out << endl;
        }
        out.close();
    }

    void save_NSadj(int process_number) {
        ofstream out(string("NSadj") + to_string(process_number) + ".msh");
        tIdx<IDX_N> ns_beg, ns_end;
        ns_beg = 0;
        ns_end = NSadj.elem_number();
        for (tIdx<IDX_N> i = ns_beg; i < ns_end; ++i) {
            out << NSadj.GetBlockSize(i);
            for (int j = 0; j < NSadj.GetBlockSize(i); ++j) {
                out << " " << NSadj[i][j].idx();
            }
            out << endl;
        }
        out.close();
    }

    void save_SNadj(int process_number) {
        ofstream out(string("SNadj") + to_string(process_number) + ".msh");
        for (tInt i = 0; i < (tInt) edges.size(); ++i) {
            out << edges[i].first << " " << edges[i].second << endl;
        }
        out.close();
    }

/// prints

    void print_file(const string& name) {
        print_file(name, name);
    }

    void print_file(const string& name, const string& title) {
        ifstream in;
        in.open(name);
        string s;
        cout << title << endl;
        while (getline(in, s)) {
            cout << "    " << s << endl;
        }
        in.close();
    }

    void save_all(int process_number) {
        save_coord(process_number);
        save_back_topo(process_number);
        save_mesh(process_number);
        save_nbf(process_number);
        save_topo(process_number);
        save_ESadj(process_number);
        save_SEadj(process_number);
        save_NSadj(process_number);
        save_SNadj(process_number);
        save_NNadj(process_number);
    }

    void gen_all() {
        double start_time = omp_get_wtime();
        gen_inverse_topo();
        cout << "generating inverse_topo " << omp_get_wtime() - start_time << endl;

        start_time = omp_get_wtime();
        gen_NNadj();
        cout << "generating NNadj " << omp_get_wtime() - start_time << endl;

        start_time = omp_get_wtime();
        gen_edges();
        cout << "generating edges " << omp_get_wtime() - start_time << endl;

        start_time = omp_get_wtime();
        gen_NSadj();
        cout << "generating NSadj " << omp_get_wtime() - start_time << endl;


        start_time = omp_get_wtime();
        gen_ESadj();
        cout << "generating ESadj " << omp_get_wtime() - start_time << endl;

        start_time = omp_get_wtime();
        gen_SEadj();
        cout << "generating SEadj " << omp_get_wtime() - start_time << endl;
    }
};

int meshGen(Mesh& m, string filePath) {
    tInt Nx, Ny, K, M;
    double Lx, Ly;
    ifstream fin;
    try {
        fin.open(filePath);
    }
    catch (...) {
        cout << "No such file" << endl;
        return 1;
    }
    fin >> Nx >> Ny >> Lx >> Ly >> K >> M;
    m.child = new Mesh;
    m.child->father = &m;

    //Nodes
    vector< tIdx<IDX_N> > allChildsNodes;
    m.nodes.clear();
    double xSize = Lx / (Nx - 1), ySize = Ly / (Ny - 1);
    for (tInt i = 0; i < Ny; ++i) {
        for (tInt j = 0; j < Nx; ++j) {
            m.nodes.push_back(Node(j * xSize, i * ySize));
            if (i == 0 || i == Ny - 1|| j == 0 || j == Nx - 1) {
                tIdx<IDX_N> childsNodes;
                childsNodes = m.child->nodes.size();
                int fathersNodes = (int) m.nodes.size() - 1;
                m.child->nodes.push_back(Node(j * xSize, i * ySize));
                allChildsNodes.push_back(childsNodes);
                m.nodeConnectionWithChild[fathersNodes] = childsNodes.idx();
                m.child->nodeConnectionWithFather.push_back(fathersNodes);
            }
        }
    }
    //Elements
    m.T.clear();
    m.child->T.clear();
    m.child->T.add(allChildsNodes);
    tInt now = 1;
    bool split = true, direction = true;
    for (tInt i = 0; i < Ny - 1; ++i) {
        for (tInt j = 0; j < Nx - 1; ++j) {
            vector< tIdx<IDX_N> > tm;
            tIdx<IDX_N> leftUp;
            leftUp = i * Nx + j;
            tIdx<IDX_N> leftDown;
            leftDown = leftUp + Nx;
            if (split) {
                if (direction) {
                    tm = {leftUp, leftDown, leftDown + 1};
                    m.T.add(tm);
                    tm = {leftUp, leftUp + 1, leftDown + 1};
                    m.T.add(tm);
                }
                else {
                    tm = {leftUp, leftUp + 1, leftDown};
                    m.T.add(tm);
                    tm = {leftUp + 1, leftDown, leftDown + 1};
                    m.T.add(tm);
                }
            }
            else {
                tm = {leftUp, leftUp + 1, leftDown, leftDown + 1};
                m.T.add(tm);
            }
            if (split && now % K == 0) {
                direction ^= true;
                now = 0;
                split = false;
            }
            else if (!split && now % M == 0) {
                now = 0;
                split = true;
            }
            ++now;
        }
    }
    //Boundary
    m.boundary.clear();
    for (tInt i = 0; i < Ny - 1; ++i) {
        Boundary b1 = Boundary(i * Nx, (i + 1) * Nx, 1), b2 = Boundary((i + 1) * Nx - 1, (i + 2) * Nx - 1, 2);
        m.boundary.push_back(b1);
        m.boundary.push_back(b2);
        m.child->boundary.push_back(b1);
        m.child->boundary.push_back(b2);
    }
    for (tInt i = 0; i < Nx - 1; ++i) {
        Boundary b1 = Boundary(i, i + 1, 3), b2 = Boundary(i + (Ny - 1) * Nx, i + 1 + (Ny - 1) * Nx, 4);
        m.boundary.push_back(b1);
        m.boundary.push_back(b2);
        m.child->boundary.push_back(b1);
        m.child->boundary.push_back(b2);
    }
    return 0;
}

void test_errors(const Mesh* mesh) {
    cout << endl << "Testing some errors:" << endl;
    tIdx<IDX_E> i;
    auto T = mesh->T;
    i = 90;
    try { T[i]; } catch (...) {}
    i = 0;
    try { T[i][3]; } catch (...) {}
    try { T[i][-1]; } catch (...) {}
}

void help() {
    cout << "Enter path to file, that contains 6 numbers:" << endl;
    cout << "Number of nodes by X" << endl << "Number of nodes by Y" << endl;
    cout << "Width" << endl << "Height" << endl;
    cout << "Number of divided cells" << endl << "Number of undivided cells" << endl;
}

int main(int argc, char* argv[]) {
    string s;
    if (argc == 1) {
        s = "in.txt";
        help();
        cout << endl << "Trying to find needed info in 'in.txt'" << endl << endl;
    } else {
        bool name = false;
        for (int i = 1; i < argc; ++i) {
            if (argv[i][0] != '-') {
                name = true;
                s = argv[i];
            }
        }
        if (!name) {
            help();
            return 0;
        }
    }
    Mesh* mesh = new Mesh;

    double start_time = omp_get_wtime();
    if (meshGen(*mesh, s)) {
        return 1;
    }
    cout << "meshGen (nodes and Topo) time " << omp_get_wtime() - start_time << endl;

    mesh->gen_all();
    cout << "all generating time " << omp_get_wtime() - start_time << endl;
    size_t mem = mesh->GetOccupiedMemory();
    cout << "Occuppied memory is about " << mem << " bytes" << endl;

    int process_number = 0;
    mesh->save_all(process_number);

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-mesh") {
            mesh->print_file("mesh" + to_string(process_number) + ".txt");
            mesh->print_file("coordinate" + to_string(process_number) + ".msh");
            mesh->print_file("topo" + to_string(process_number) + ".msh");
            mesh->print_file("bctopo" + to_string(process_number) + ".msh");
        } else if (arg == "-EN") {
            mesh->print_file("topo" + to_string(process_number) + ".msh", "EN");
        } else if (arg == "-NE") {
            mesh->print_file("backtopo" + to_string(process_number) + ".msh", "NE");
        } else if (arg == "-NN") {
            mesh->print_file("coordinate" + to_string(process_number) + ".msh", "NN");
        } else if (arg == "-SN") {
            mesh->print_file("SNadj" + to_string(process_number) + ".msh", "SN");
        } else if (arg == "-NS") {
            mesh->print_file("NSadj" + to_string(process_number) + ".msh", "NS");
        } else if (arg == "-SE") {
            mesh->print_file("SEadj" + to_string(process_number) + ".msh", "SE");
        } else if (arg == "-ES") {
            mesh->print_file("ESadj" + to_string(process_number) + ".msh", "ES");
        }
    }

    bool debug_output = false;
    if (debug_output) {
        mesh->print_all_debug();
    }

    bool test_output = false;
    if (test_output) {
        test_errors(mesh);
    }

    cout << "Done" << endl;
    return 0;
}
