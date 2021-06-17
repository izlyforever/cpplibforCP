<head>
	<script type="text/x-mathjax-config">
		MathJax.Hub.Config({
		  tex2jax: {
			skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
			inlineMath: [['$','$']],
			processEscapes: true
		  }
		});
	</script>
	<script type="text/javascript" async
	  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-MML-AM_CHTML">
	</script>
</head>

# cpplibForCP

C++17 (-O2) template for competitive programming algorithms, which contains numerous math algorithms.

- Aims: build a stable, fast, easy-to-read C++ template.(the more you use it, the more reliable it is)
- Document: [English version](document_en) and [Chinese version](document_cn) for addition.
- Reference: [CP-algorithm](https://cp-algorithms.com/), [AtCoder](https://github.com/atcoder/ac-library) and many awesome blogs, for instance [zscoder](https://codeforces.com/profile/zscoder) posted in codeforces.
- Additional Chinese Reference: [OI-wiki](https://oi-wiki.org/) and [LuoGu](https://www.luogu.com.cn/)
- Test: [Welcome Your Test](test)

> It's really hard to give a reference for every algorithm (in English).


## Need help

I don't know how to enable mathjax in github Page without adding origin source of jekyll. 

so I add following code in every `.md` file

``` html
<head>
	<script type="text/x-mathjax-config">
		MathJax.Hub.Config({
		  tex2jax: {
			skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
			inlineMath: [['$','$']],
			processEscapes: true
		  }
		});
	</script>
	<script type="text/javascript" async
	  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-MML-AM_CHTML">
	</script>
</head>
```