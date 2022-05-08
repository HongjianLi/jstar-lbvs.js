import path from 'path';
import util from 'util';
import { read, promises as fs } from 'fs';
import local_time_string from './utility.js';

/**
 * Represents a compound database.
 */
class CompoundDatabase {
	/**
	 * @returns {Array<string>} Descriptor files.
	 */
	static get descriptors() {
		return ['natm.u16', 'nhbd.u16', 'nhba.u16', 'nrtb.u16', 'nrng.u16', 'xmwt.f32', 'tpsa.f32', 'clgp.f32', 'usrcat.f32', 'conformers.sdf.ftr.u64'];
	}

	/**
	 * @param {string} dpth - Path to the compound database directory.
	 */
	constructor(dpth) {
		this.dpth = dpth;
		this.name = path.basename(dpth);
	}

	/**
	 * Read descriptor files.
	 */
	async read_descriptors() {
		console.log(`${local_time_string()} Reading ${this.name}`);

		// Read molecular descriptor files.
		const { descriptors } = CompoundDatabase;
		for (let k = 0; k < descriptors.length; ++k) {
			const filename = descriptors[k];
			const buf = await fs.readFile(path.join(this.dpth, filename));
			console.log(`${local_time_string()} Read ${buf.length} bytes from ${filename}`);
			switch (filename.slice(-3)) {
				case 'u16':
					this[filename] = new Uint16Array(buf.buffer, 0, buf.length / Uint16Array.BYTES_PER_ELEMENT);
					break;
				case 'f32':
					this[filename] = new Float32Array(buf.buffer, 0, buf.length / Float32Array.BYTES_PER_ELEMENT);
					break;
				case 'u64':
					this[filename] = new BigUint64Array(buf.buffer, 0, buf.length / BigUint64Array.BYTES_PER_ELEMENT);
					break;
				default:
					console.error(`Unsupported descriptor file extension: ${filename}`);
			}
		}
		this.num_compounds = this[descriptors[0]].length;
		descriptors.slice(1, -2).forEach(filename => { // except 'usrcat.f32' and 'conformers.sdf.ftr.u64'
			console.assert(this[filename].length === this.num_compounds);
		});
		this.num_conformers = this['conformers.sdf.ftr.u64'].length;
		console.assert(this.num_conformers === this['usrcat.f32'].length / 60);
		console.assert(this.num_conformers === this.num_compounds << 2);
		console.log(`${local_time_string()} Found ${this.num_compounds} compounds and ${this.num_conformers} conformers`);

		// Open conformers.sdf to obtain a file descriptor.
		this['conformers.sdf'] = await fs.open(path.join(this.dpth, "conformers.sdf"));
	}

	/**
	 * Read the ith conformer out of conformers.sdf
	 * @param {number} index - Index of the conformer to read.
	 * @returns {string} Conformer string.
	 */
	async read_conformer(index) {
		const position = index ? this['conformers.sdf.ftr.u64'][index - 1] : 0n;
		const length = parseInt(this['conformers.sdf.ftr.u64'][index] - position);
		const buffer = Buffer.alloc(length);
		const { bytesRead } = await readPromisied(this['conformers.sdf'].fd, { buffer, position }); // fs.read(fd[, options], callback) supports position <bigint>. https://nodejs.org/api/fs.html#fsreadfd-options-callback
		console.assert(bytesRead === length);
		return buffer.toString();
	}
}

const readPromisied = util.promisify(read); // Convert the callback style of fs.read to the promise style.

export { CompoundDatabase };
