# cpplibForCP

C++17 (-O2) template for competitive programming algorithms, which contains numerous math algorithms.

- Aims: build a stable, fast, easy-to-read C++ template.(the more you use it, the more reliable it is)
- Document: [English version](https://izlyforever.github.io/cpplibforCP/index) and many comments in the code, [Chinese version](https://izlyforever.github.io/cpplibforCP/cn) for addition.
- Reference: [CP-algorithm](https://cp-algorithms.com/), [AtCoder](https://github.com/atcoder/ac-library) and many awesome blogs, for instance [zscoder](https://codeforces.com/profile/zscoder) posted in codeforces and [jiangly's submissions](https://codeforces.com/submissions/jiangly)
- Additional Chinese Reference: [OI-wiki](https://oi-wiki.org/) and [LuoGu](https://www.luogu.com.cn/)
- Test: almost all algorithms are tested in [LuoGu](https://www.luogu.com.cn/) and [Codeforces](https://codeforces.com/)
- Mkdocs: The Document page is build with MkDocs

> It's really hard to for me to give a reference for every algorithm (in English).

## ToDo

- Add GTest
- support CI/CD
- Refine factor using `last spf`


## Deploy

``` shell
mkdocs gh-deploy
cloudbase hosting:deploy site cpplibforCP -e  blog-5g9vf63s2e355403
```

### Requirement(may need venv)

- `pip install mkdocs`
- `pip install mkdocs-bootswatch`
- `pip install python-markdown-math`
