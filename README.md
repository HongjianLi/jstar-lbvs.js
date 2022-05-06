# jstar-lbvs.js
[jstar]'s daemon for ligand-based virtual screening (LBVS), written in JavaScript.

A C++ implementation is available at https://github.com/HongjianLi/jstar-lbvs.cpp

## Usage
    node main.js --databases ../jstar/databases --host localhost --port 27017 --user jstard --pass password

## Architecture
![jstar architecture](https://github.com/HongjianLi/jstar/blob/master/public/architecture.png)

## Components
### Database
* [MongoDB]
### Daemon
* [MongoDB Node.js Driver]
* [rdkit]
* [piscina]
* [commander]

## License
[MIT License]

## Developers
* [Jacky Lee]
* Jessie Sze

## Logo
![jstar logo](https://github.com/HongjianLi/jstar/blob/master/public/logo.svg)

[jstar]: https://github.com/HongjianLi/jstar
[MongoDB]: https://github.com/mongodb/mongo
[MongoDB Node.js Driver]: https://mongodb.github.io/node-mongodb-native
[rdkit]: https://github.com/rdkit/rdkit/tree/master/Code/MinimalLib
[piscina]: https://github.com/piscinajs/piscina
[commander]: https://github.com/tj/commander.js
[MIT License]: https://github.com/HongjianLi/jstar-lbvs.js/blob/master/LICENSE
[Jacky Lee]: https://github.com/HongjianLi
