## Metagenomix

Metagenomix is a Python module for analyzing species composition of a sequenced sample of DNA or RNA.

## Code Example


Show what the library does as concisely as possible, developers should be able to figure out **how** your project solves their problem by looking at the code example. Make sure the API you are showing off is obvious, and that your code is short and concise.

## Motivation

A short description of the motivation behind the creation and maintenance of the project. This should explain **why** the project exists.

## Installation

To install Metagenomix, either clone or download this repository and run:
```
python setup.py build install
```
After you've built Metagenomix, you need to download the information on taxa names, ranks and relationships:
```
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
dmp2taxtree nodes.dmp names.dmp ./metagenomix/
python setup.py build install
```
After the build step, you need to setup the gi to taxonomy mapping database which will enable Metagenomix to map particular organism names to your GI identifiers.


## Contributors

This project is being developed my Ana Bulović.
If you have any questions, or wish to leave a comment, you can contact me on https://twitter.com/praxaton, or using a contact form on http://metagen.zesoi.fer.hr.

## License

Copyright (c) 2015 Ana Bulović metagen.zesoi.fer.hr

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
