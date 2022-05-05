import path from 'path';
import fs from 'fs';

/**
 * Represents a compound database.
 */
class compound_database {
	/**
	 * compound_database::compound_database(const path dpth)
	 * @param {string} dpth - Path to the compound database directory.
	 */
	constructor(dpth) {
		/**
		 * Path to the compound database directory.
		 * @type {string}
		 */
		this.dpth = dpth;

		/**
		 * Database name.
		 * @type {string}
		 */
		this.name = path.basename(dpth);
		console.log(`Reading ${this.name}`);

		// Read molecular descriptor files.
		/**
		 * Number of atoms.
		 * vector<uint16_t> natm;
		 * @type {Array<number>}
		 */
		const natm = fs.readFileSync(path.join(dpth, "natm.u16"));
		this.natm = new Uint16Array(natm.buffer, 0, natm.length / Uint16Array.BYTES_PER_ELEMENT);
		this.num_compounds = this.natm.length;
		/**
		 * Number of hydrogen bond donors.
		 * vector<uint16_t> nhbd;
		 * @type {Array<number>}
		 */
		const nhbd = fs.readFileSync(path.join(dpth, "nhbd.u16"));
		this.nhbd = new Uint16Array(nhbd.buffer, 0, nhbd.length / Uint16Array.BYTES_PER_ELEMENT);
		console.assert(this.nhbd.length === this.num_compounds);
		/**
		 * Number of hydrogen bond acceptors.
		 * vector<uint16_t> nhba;
		 * @type {Array<number>}
		 */
		const nhba = fs.readFileSync(path.join(dpth, "nhba.u16"));
		this.nhba = new Uint16Array(nhba.buffer, 0, nhba.length / Uint16Array.BYTES_PER_ELEMENT);
		console.assert(this.nhba.length === this.num_compounds);
		/**
		 * Number of rotatable bonds.
		 * vector<uint16_t> nrtb;
		 * @type {Array<number>}
		 */
		const nrtb = fs.readFileSync(path.join(dpth, "nrtb.u16"));
		this.nrtb = new Uint16Array(nrtb.buffer, 0, nrtb.length / Uint16Array.BYTES_PER_ELEMENT);
		console.assert(this.nrtb.length === this.num_compounds);
		/**
		 * Number of rings.
		 * vector<uint16_t> nrng;
		 * @type {Array<number>}
		 */
		const nrng = fs.readFileSync(path.join(dpth, "nrng.u16"));
		this.nrng = new Uint16Array(nrng.buffer, 0, nrng.length / Uint16Array.BYTES_PER_ELEMENT);
		console.assert(this.nrng.length === this.num_compounds);
		/**
		 * Exact molecular weight.
		 * vector<float> xmwt;
		 * @type {Array<number>}
		 */
		const xmwt = fs.readFileSync(path.join(dpth, "xmwt.f32"));
		this.xmwt = new Float32Array(xmwt.buffer, 0, xmwt.length / Float32Array.BYTES_PER_ELEMENT);
		console.assert(this.xmwt.length === this.num_compounds);
		/**
		 * Topological polar surface area.
		 * vector<float> tpsa;
		 * @type {Array<number>}
		 */
		const tpsa = fs.readFileSync(path.join(dpth, "tpsa.f32"));
		this.tpsa = new Float32Array(tpsa.buffer, 0, tpsa.length / Float32Array.BYTES_PER_ELEMENT);
		console.assert(this.tpsa.length === this.num_compounds);
		/**
		 * clogP
		 * vector<float> clgp;
		 * @type {Array<number>}
		 */
		const clgp = fs.readFileSync(path.join(dpth, "clgp.f32"));
		this.clgp = new Float32Array(clgp.buffer, 0, clgp.length / Float32Array.BYTES_PER_ELEMENT);
		console.assert(this.clgp.length === this.num_compounds);

		// Read usrcat feature file.
		/**
		 * USRCAT features.
		 * vector<array<float, 60>> usrcat;
		 * @type {Array<number>}
		 */
		const usrcat = fs.readFileSync(path.join(dpth, "usrcat.f32"));
		this.usrcat = new Float32Array(usrcat.buffer, 0, usrcat.length / Float32Array.BYTES_PER_ELEMENT);
		this.num_conformers = this.usrcat.length / 60;
		console.assert(this.num_conformers === this.num_compounds << 2);

		// Read conformers.sdf.ftr footer file.
		/**
		 * Footer file of conformers.sdf
		 * vector<size_t> conformers_sdf_ftr;
		 * @type {Array<number>}
		 */
		const conformers_sdf_ftr = fs.readFileSync(path.join(dpth, "conformers.sdf.ftr"));
		this.conformers_sdf_ftr = new BigUint64Array(conformers_sdf_ftr.buffer, 0, conformers_sdf_ftr.length / BigUint64Array.BYTES_PER_ELEMENT);
		console.assert(this.conformers_sdf_ftr.length === this.num_conformers);

		// Open conformers.sdf to obtain a file descriptor.
		this.conformers_sdf = fs.openSync(path.join(this.dpth, "conformers.sdf"));
		console.log(`Found ${this.num_compounds} compounds and ${this.num_conformers} conformers`);
	}

	/**
	 * Read the ith conformer out of conformers.sdf
	 * string compound_database::read_conformer(const size_t index, ifstream& ifs) const
	 * @param {*} index 
	 */
	read_conformer(index) {
		const position = index ? this.conformers_sdf_ftr[index - 1] : 0n;
		const length = parseInt(this.conformers_sdf_ftr[index] - position);
		const buffer = Buffer.alloc(length);
		const bytesRead = fs.readSync(this.conformers_sdf, buffer, { position });
		console.assert(bytesRead === length);
		return buffer.toString();
	}
}

export { compound_database };