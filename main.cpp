#include<bits/stdc++.h>
using namespace std;

struct Data{
    int numItems;
    vector<double> dataItems;

    Data(int num, vector<double> items){
        this->numItems = num;
        for(int i=0;i<num;i++){
            this->dataItems.push_back(items[num]);
        }
    } 
};

class Compare{
public:
    bool operator() (q_elem a, q_elem b){
        if(a.second <= b.second){
            return true;
        }else if(a.second == b.second){
            return a.first <= b.first;
        }else{
            return false;
        }
    }
};

typedef pair<int, double> q_elem;
typedef priority_queue<q_elem, vector<q_elem>, Compare> dataQueue;

double cosine_dist(Data a, Data b){
    assert(a.numItems == b.numItems);
    double dotProd = 0.0;
    double lenA = 0.0;
    double lenB = 0.0;
    for(int i=0;i<a.numItems;i++){
        double ad = a.dataItems[i];
        double bd = b.dataItems[i];
        dotProd += (ad * bd);
        lenA += (ad * ad);
        lenB += (bd * bd);
    }
    dotProd = 1.0 - dotProd / (sqrt(lenA) * sqrt(lenB));
    return dotProd;
}

void trim(dataQueue &dq, int limit){
    vector<q_elem> v;
    while(dq.size()>0){
        v.push_back(dq.top());
        dq.pop();
    }
    for(int i=0;i<limit;i++){
        dq.push(v[i]);
    }
}

void SearchLayer(Data q, vector<int> indptr, vector<int> index, vector<int> level_offset, int currLevel, unordered_set<int> &visited, vector<Data> vect, dataQueue &candidates){
    dataQueue *newCandidates = new dataQueue(candidates);
    int k = candidates.size();
    while(newCandidates->size() > 0){
        int curr = newCandidates->top().first;
        newCandidates->pop();
        int start = indptr[curr] + level_offset[currLevel];
        int end = indptr[curr] + level_offset[currLevel + 1];
        for(int i = start;i<end;i++){
            int currIdx = index[i];
            if(visited.find(currIdx) == visited.end() || currIdx == -1){
                continue;
            }
            visited.insert(currIdx);
            double currDist = cosine_dist(q, vect[currIdx]);
            //Do distance check
            candidates.push(make_pair(currIdx, currDist));
            trim(candidates,k);
            newCandidates->push(make_pair(currIdx,currDist));
        }
    }
}

dataQueue queryHNSW(Data q, int top_k, int ep, vector<int> indptr, vector<int> index, vector<int> level_offset, int max_level, vector<Data> vect){
    dataQueue candidates;
    candidates.push(make_pair(ep, cosine_dist(q, vect[ep])));
    unordered_set<int> visited;
    visited.insert(ep);
    for(int i = max_level - 1;i>=0;i--){
        SearchLayer(q, indptr, index, level_offset, i, visited, vect, candidates);
    } 
    return candidates;
}

int main(){

}