import path from 'path';
import { readSync, promises as fs } from 'fs';
import local_time_string from './utility.js';

/**
 * Represents a compound database.
 */
class compound_database {
	static get descriptors() {
		return {
			u16: ['natm', 'nhbd', 'nhba', 'nrtb', 'nrng'],
			f32: ['xmwt', 'tpsa', 'clgp', 'usrcat'],
			u64: ['conformers.sdf.ftr'],
		};
	}
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
	}

	async read_descriptors() {
		console.log(`${local_time_string()} Reading ${this.name}`);

		// Read molecular descriptor files.
		const { u16, f32, u64 } = compound_database.descriptors;
		for (let k = 0; k < u16.length; ++k) {
			const key = u16[k];
			const filename = `${key}.u16`;
			const buf = await fs.readFile(path.join(this.dpth, filename));
			console.log(`${local_time_string()} Read ${buf.length} bytes from ${filename}`);
			this[key] = new Uint16Array(buf.buffer, 0, buf.length / Uint16Array.BYTES_PER_ELEMENT);
		}
		for (let k = 0; k < f32.length; ++k) {
			const key = f32[k];
			const filename = `${key}.f32`;
			const buf = await fs.readFile(path.join(this.dpth, filename));
			console.log(`${local_time_string()} Read ${buf.length} bytes from ${filename}`);
			this[key] = new Float32Array(buf.buffer, 0, buf.length / Float32Array.BYTES_PER_ELEMENT);
		}
		for (let k = 0; k < u64.length; ++k) {
			const key = u64[k];
			const filename = `${key}.u64`;
			const buf = await fs.readFile(path.join(this.dpth, filename));
			console.log(`${local_time_string()} Read ${buf.length} bytes from ${filename}`);
			this[key] = new BigUint64Array(buf.buffer, 0, buf.length / BigUint64Array.BYTES_PER_ELEMENT);
		}
		this.num_compounds = this[u16[0]].length;
		u16.slice(1).concat(f32.slice(0, -1)).forEach(key => {
			console.assert(this[key].length === this.num_compounds);
		});
		this.num_conformers = this['conformers.sdf.ftr'].length;
		console.assert(this.num_conformers === this.usrcat.length / 60);
		console.assert(this.num_conformers === this.num_compounds << 2);
		console.log(`${local_time_string()} Found ${this.num_compounds} compounds and ${this.num_conformers} conformers`);

		// Open conformers.sdf to obtain a file descriptor.
		this['conformers.sdf'] = await fs.open(path.join(this.dpth, "conformers.sdf"));
	}

	/**
	 * Read the ith conformer out of conformers.sdf
	 * string compound_database::read_conformer(const size_t index, ifstream& ifs) const
	 * @param {*} index 
	 */
	read_conformer(index) {
		const position = index ? this['conformers.sdf.ftr'][index - 1] : 0n;
		const length = parseInt(this['conformers.sdf.ftr'][index] - position);
		const buffer = Buffer.alloc(length);
		const bytesRead = readSync(this['conformers.sdf'].fd, buffer, { position }); // fs.readSync(fd, buffer[, options]) supports position <bigint>.
		console.assert(bytesRead === length);
		return buffer.toString();
	}
}

export { compound_database };