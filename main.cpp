#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
using namespace std;

using q_elem = pair<int, float>;

struct Compare{
    bool operator() (q_elem const a, q_elem const b) {
        if(a.second < b.second) {
            return true;
        }
        return (a.second == b.second) && a.first < b.first;
    }
};

using dataQueue = priority_queue<q_elem, vector<q_elem>, Compare>;

float cosine_dist(vector<float> a, float *b){
    int n = a.size();
    float dotProd = 0.0;
    float lenA = 0.0;
    float lenB = 0.0;
    for(int i = 0; i < n; i++) {
        float ad = a[i];
        float bd = b[i];
        dotProd += (ad * bd);
        lenA += (ad * ad);
        lenB += (bd * bd);
    }
    dotProd = 1.0 - dotProd / (sqrt(lenA) * sqrt(lenB));
    return dotProd;
}

void trim(dataQueue &dq, int limit){
    while (dq.size() > limit)
        dq.pop();
}

dataQueue maxToMinOrViceVersaHeap(dataQueue dq) {
    vector<q_elem> v;
    dataQueue newDq;
    while (!dq.empty()) {
        q_elem top = dq.top();
        newDq.push({-top.first, -top.second});
        dq.pop();
    }
    return newDq;
}

float **vect;
vector<vector<float>> users;

void SearchLayer(vector<float> &q, int k, vector<uint32_t> &indptr, vector<int32_t> &index, vector<uint32_t> &level_offset, int currLevel, unordered_set<int> &visited, dataQueue &candidates) {
    dataQueue newCandidates = maxToMinOrViceVersaHeap(candidates);
    while(newCandidates.size() > 0) {
        int curr = -newCandidates.top().first;
        newCandidates.pop();
        int start = indptr[curr] + level_offset[currLevel];
        int end = indptr[curr] + level_offset[currLevel + 1];
        for(int i = start;i<end;i++){
            int currIdx = index[i];
            if(visited.count(currIdx) || currIdx == -1) {
                continue;
            }
            visited.insert(currIdx);
            float currDist = cosine_dist(q, vect[currIdx]);
            candidates.push(make_pair(currIdx, currDist));
            trim(candidates, k);
            newCandidates.push(make_pair(-currIdx, -currDist));
        }
    }
}

dataQueue queryHNSW(vector<float> &q, int top_k, int ep, vector<uint32_t> &indptr, vector<int32_t> &index, vector<uint32_t> &level_offset, int max_level){
    dataQueue candidates;
    candidates.push(make_pair(ep, cosine_dist(q, vect[ep])));
    unordered_set<int> visited;
    visited.insert(ep);
    for(int i = max_level;i>=0;i--) {
        SearchLayer(q, top_k, indptr, index, level_offset, i, visited, candidates);
    }
    return candidates;
}

ifstream open_file(string filename){
    ifstream file;
    file.open(filename, ios::binary);
    if(!file.is_open()){
        cout << "Could not open file" << endl;
        exit(1);
    }
    return file;
}

template <typename T, size_t S>
void read_val(T *n, ifstream &file) {
    unsigned char buffer[S];
    file.read((char *)buffer, S);
    *n = *(T *)buffer;
}

int main(int argc, char *argv[]){
    assert(argc == 5);

    string data_dir = argv[1];
    int top_k = stoi(argv[2]);
    string user = argv[3];
    string out_file = argv[4];

    int rank, size;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ifstream params_file = open_file(data_dir + "/params");
    uint32_t ep, max_level, l, d;
    read_val<uint32_t, 4>(&ep, params_file);
    read_val<uint32_t, 4>(&max_level, params_file);
    read_val<uint32_t, 4>(&l, params_file);
    read_val<uint32_t, 4>(&d, params_file);
    params_file.close();
    cout << "ep: " << ep << ' ' << "max_level: " << max_level << ' ' << "l: " << l << ' ' << "d: " << d << endl;

    ifstream indptr_file = open_file(data_dir + "/indptr");
    vector<uint32_t> indptr;
    while (indptr_file.peek() != EOF) {
        uint32_t curr;
        read_val<uint32_t, 4>(&curr, indptr_file);
        indptr.push_back(curr);
    }
    indptr_file.close();
    cout << "indptr read" << endl;

    ifstream index_file = open_file(data_dir + "/index");
    vector<int32_t> index;
    while (index_file.peek() != EOF) {
        int32_t curr;
        read_val<int32_t, 4>(&curr, index_file);
        index.push_back(curr);
    }
    index_file.close();
    cout << "index read" << endl;

    ifstream level_offset_file = open_file(data_dir + "/level_offset");
    vector<uint32_t> level_offset;
    while (level_offset_file.peek() != EOF) {
        uint32_t n;
        read_val<uint32_t, 4>(&n, level_offset_file);
        level_offset.push_back(n);
    }
    level_offset_file.close();
    cout << "level_offset read" << endl;

    ifstream vect_file = open_file(data_dir + "/vect");
    vect = new float*[l];
    for(int i = 0; i < l; i++) {
        vect[i] = new float[d];
    }
    for (int j = 0; j < l; j++) {
        for (int i = 0; i < d; i++) {
            read_val<float, 4>(&vect[j][i], vect_file);
        }
    }
    vect_file.close();
    cout << "vect read" << endl;

    ifstream user_file(user);
    if (!user_file.is_open()) {
        cout << "Could not open file" << endl;
        exit(1);
    }
    while (user_file.peek() != EOF) {
        vector<float> items;
        for (int i = 0; i < d; i++) {
            if (user_file.peek() == EOF) {
                break;
            }
            float dbl;
            user_file >> dbl;
            items.push_back(dbl);
        }
        if (items.size() < d) {
            break;
        }
        users.push_back(items);
    }
    user_file.close();
    cout << "users read" << endl;

    int num_threads[size];
    int user_start_indices[size + 1], user_sizes[size];
    num_threads[rank] = omp_get_max_threads();
    MPI_Allgather(&num_threads[rank], 1, MPI_INT, num_threads, 1, MPI_INT, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < size; i++)
            cout << num_threads[i] << ' ';
        cout << endl;
    }

    int num_threads_total = 0;
    for(int i = 0;i<size;i++){
        num_threads_total += num_threads[i];
    }

    int num_users = users.size();
    int num_blocks = num_users / num_threads_total;
    int remaining = num_users % num_threads_total;
    for (int node = 0; node < size; node++) {
        user_sizes[node] = num_blocks * num_threads[node];
        int extra = (remaining + size - 1) / size;
        if(node == size - 1){
            user_sizes[node] += (remaining - extra * (size - 1));
        }else{
            user_sizes[node] += extra; 
        }
        user_sizes[node] *= sizeof(q_elem) * top_k;
    }

    user_start_indices[0] = 0;
    for (int node = 1; node <= size; node++) {
        user_start_indices[node] = user_start_indices[node-1] + user_sizes[node-1];
    }
    q_elem output[num_users][top_k];

    int startNode = user_start_indices[rank] / (sizeof(q_elem) * top_k);
    int endNode = user_start_indices[rank+1] / (sizeof(q_elem) * top_k);
    #pragma omp parallel for shared(users, top_k, ep, indptr, index, level_offset, max_level, l, d, vect, output, startNode, endNode)
        for (int i=startNode; i < endNode; i++) {
            dataQueue candidates = queryHNSW(users[i], top_k, ep, indptr, index, level_offset, max_level);
            int j = 0;
            while (!candidates.empty()) {
                output[i][j++] = make_pair(candidates.top().first,candidates.top().second);
                candidates.pop();
            }
            reverse(output[i], output[i] + j);
            while (j < top_k)
                output[i][j++] = make_pair(-1, -1);
        }

    MPI_Gatherv(output[user_start_indices[rank] / (sizeof(q_elem) * top_k)], user_sizes[rank], MPI_BYTE, output, user_sizes, user_start_indices, MPI_BYTE, 0, MPI_COMM_WORLD);


    if(rank == 0) {
        ofstream out(out_file);
        for(int i=0;i<users.size();i++){
            for(int j=0;j<top_k;j++){
                if (output[i][j].first == -1)
                    break;
                out << output[i][j].first << " ";
            }
            out << '\n';
        }
        out.close();
    }

    MPI_Finalize();
    
    return 0;
}