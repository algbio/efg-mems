 /* Wrapper using part of br-index API 
  * with the bdbwt as needed for MEM finding
  */

#include "../br-index/src/definitions.hpp"
#include "../bdbwt/include/BD_BWT_index.hh" 

namespace bri {

// sample maintained during the search
struct bd_sample {
    /*
     * state variables for left_extension & right_extension
     * range: SA range of P
     * rangeR: correspondents to range in reversed text
     */
    range_t range, rangeR;
    
    bd_sample(): range(), rangeR() {}

    bd_sample(range_t range_, 
              range_t rangeR_) 
              :
              range(range_),
              rangeR(rangeR_) {}

    void set_values(range_t range_, 
                    range_t rangeR_)
    {
        range = range_;
        rangeR = rangeR_;
    }

    // the pattern does not exist
    bool is_invalid()
    {
        return (range.first > range.second) || (rangeR.first > rangeR.second);
    }

    // range size
    ulint size()
    {
        return range.second + 1 - range.first;
    }
};

class bdbwt_index {

public:


    bdbwt_index() {}

    bdbwt_index(std::string const& input)
    {
       bdbwt = BD_BWT_index<>((uint8_t*)input.c_str());
       saidx = new ulint[bdbwt.size()];
       ulint j = LF(0);
       for (ulint i=0; i<bdbwt.size()-2; i++) {
          saidx[j] = bdbwt.size()-i-2; 
          j = LF(j);
       }
    }
    
     /*
     * get full BWT range
     */
    range_t full_range()
    {
        return {0,bdbwt.size()-1};
    }


    /*
     * get a sample corresponding to an empty string
     */
    bd_sample get_initial_sample(bool reverse = false)
    {
       return bd_sample(full_range(), // entire SA range
                            full_range()); // entire SAR range   
    }
     
    ulint LF(ulint j)
    {
       return bdbwt.backward_step(j);
    }

    ulint FL(ulint j)
    {
       return bdbwt.forward_step(j);
    }
    
      
    /*
     * search the pattern cP (P:the current pattern)
     * returns SA&SAR range corresponding to cP
     */
    bd_sample left_extension(uchar c, bd_sample const& prev_sample)
    {
        Interval_pair pair = bdbwt.left_extend(Interval_pair(prev_sample.range.first,prev_sample.range.second,prev_sample.rangeR.first,prev_sample.rangeR.second),c);
        bd_sample sample;
        sample.range.first=pair.forward.left;
        sample.range.second=pair.forward.right;
        sample.rangeR.first=pair.reverse.left;
        sample.rangeR.second=pair.reverse.right;
        return sample;
    }

    /*
     * search the pattern Pc (P:the current pattern)
     * return SAR&SA range corresponding to Pc
     */
    bd_sample right_extension(uchar c, bd_sample const& prev_sample)
    {
        Interval_pair pair = bdbwt.right_extend(Interval_pair(prev_sample.range.first,prev_sample.range.second,prev_sample.rangeR.first,prev_sample.rangeR.second),c);
        bd_sample sample;
        sample.range.first=pair.forward.left;
        sample.range.second=pair.forward.right;
        sample.rangeR.first=pair.reverse.left;
        sample.rangeR.second=pair.reverse.right;
        return sample;
    }


    /*
     * get BWT[i] or BWT^R[i]
     */
    uchar bwt_at(ulint i, bool reversed = false)
    {
        if (!reversed) bdbwt.forward_bwt_at(i);
        return bdbwt.backward_bwt_at(i);
    }



    /*
     * save index to "{path_prefix}.bdbwt" file
     */
    void save_to_file(std::string const& path_prefix)
    {

        std::string path = path_prefix + ".bdbwt";
        
        bdbwt.save_to_disk(path,"");
    
    }

    /*
     * load index file from path
     */
    void load_from_file(std::string const& path)
    {

        bdbwt.load_from_disk(path,"");

    }
    
    void load(std::istream& in) {
       // not used, defined for compatibility
    }

    ulint text_size() { return bdbwt.size() - 1; }

    ulint bwt_size(bool reversed=false) { return bdbwt.size(); }

    uchar get_terminator() {
        return bdbwt.END;
    }

    ulint get_terminator_position(bool reversed = false)
    {
	    return FL(0);
    }
    
    std::vector<ulint> locate_sample(bd_sample const& sample)
    {
       ulint n_occ = sample.range.second + 1 - sample.range.first;
       std::vector<ulint> res(n_occ);
       for (ulint i=sample.range.first; i <= sample.range.second; i++)
          res[i-sample.range.first] = saidx[i];
       return res; 
    }

    /*
     * get string representation of BWT
     */
    std::string get_bwt(bool reversed = false)
    {
        if (!reversed)
        {
            std::string res(bdbwt.size(), '0');
            for (size_t i = 0; i < bdbwt.size(); ++i)
                res[i] = bdbwt.forward_bwt_at(i);
            return res;
        } else {
            std::string res(bdbwt.size(), '0');
            for (size_t i = 0; i < bdbwt.size(); ++i)
                res[i] = bdbwt.backward_bwt_at(i);
            return res;
        }
    }
    


private:


   BD_BWT_index<> bdbwt;
   ulint* saidx;

};

};


