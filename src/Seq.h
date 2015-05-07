#ifndef SEQ_H
#define SEQ_H

class Seq {
    public:
        Seq();
        Seq(const char *seq_str, size_t len, const std::string& description);
        Seq(const Seq& seq);
        ~Seq();
        Seq& operator= (const Seq& seq);

        size_t len;
        const std::string& desc() const { return description; }

        std::string to_string() const;

        //std::string desc;
        uint8_t *data;

        uint32_t substr(size_t pos, size_t len) const;

    private:
        std::string description;
        unsigned int bytes;
};

#endif
