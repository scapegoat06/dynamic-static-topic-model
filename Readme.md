## Dynamic and Static Topic Models
We propose a dynamic and static topic model, which simultaneously considers the dynamic structures of the temporal topic evolution and the static structures of the topic hierarchyat each time. 
[paper](https://www.aclweb.org/anthology/P18-2082).

## How to run
Clone & Build:
```
git clone 〜〜〜
cd src/
g++ -Wall -O2 -o dstm.cc
g++ -Wall -O2 -c main.cc
g++ -Wall -O2 -o dstm dstm.o main.o
```

To Model Run:
```
cd src/
./dstm train.txt test.txt vocab.txt SuperTopicNumber SubTopicNumber Iterations

# example
./dstm ../data/nips_tarin.txt ../data/nips_test.txt ../data/vocab_nips.txt 5 10 50

```

- `train.txt` and  `test.txt` is formated as below.

``` 
13　//Time Step Number
90,95,101,143,144,127,144,140,152,152,151,151,150 //Each Time Step Document Size
11463 // Vocabulary Size
917702 // Total Token Size
1 1 8193 3 // time_step,doc_id,word_id,word_count
・・・・
・・・・
```

- `vocab.txt` is vocabulary of corpus

OutputFiles:

some parameter and perplexity(ppl) is written to file.

- phi:
- beta:
- alpha1:
- alpha2:
- thetas:
- thetask:
- ppl:

## Reference
If you use anything in this repository, please cite:

Rem Hida, Naoya Takeishi, Takehisa Yairi, and Koichi Hori, “Dynamic and Static Topic Model for Analyzing Time-Series Document Collections,” in Proceedings of the 56th Annual Meeting of the Association for Computational Linguistics (ACL), 2018.

```
@inproceedings{hida-etal-2018-dynamic,
    title = "Dynamic and Static Topic Model for Analyzing Time-Series Document Collections",
    author = "Hida, Rem and Takeishi, Naoya  and Yairi, Takehisa  and Hori, Koichi",
    booktitle = "Proceedings of the 56th Annual Meeting of the Association for Computational Linguistics (Volume 2: Short Papers)",
    month = jul,
    year = "2018",
    address = "Melbourne, Australia",
    publisher = "Association for Computational Linguistics",
    url = "https://www.aclweb.org/anthology/P18-2082",
    doi = "10.18653/v1/P18-2082",
    pages = "516--520"
}
```
