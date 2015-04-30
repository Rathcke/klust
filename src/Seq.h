#ifndef SEQ_H
#define SEQ_H

class Seq {
    public:
        Seq();
        Seq(const char *seq_str, size_t len, const std::string& description);
        Seq(const Seq& seq);
        ~Seq();
        Seq& operator= (const Seq& seq);

        inline uint8_t* data() const { return seq_data; }
        inline size_t length() const { return seq_len; }
        //std::string& desc() { return description; }

        std::string to_string() const;

        std::string desc;

        uint32_t substr(size_t pos, size_t len) const;

    private:
        uint8_t *seq_data;
        size_t seq_len;
        //std::string description;
        unsigned int bytes;
};

#endif
