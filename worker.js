/**
 * 
 * @param {*} param0 
 */
export function calculate({ chunk_size, l, num_compounds, scores, usrcat, qnu0, q, cnfids, zcase, num_hits }) {
	// Loop over compounds of the current chunk.
	const chunk_beg = chunk_size * l;
	const chunk_end = Math.min(chunk_beg + chunk_size, num_compounds);
	for (let k = chunk_beg; k < chunk_end; ++k) {
		// Loop over conformers of the current compound and calculate their primary score.
		for (let j = k << 2; j < ((k + 1) << 2); ++j) {
			const o = 60 * j; // offset to usrcat
			let s = 0;
			for (let i = 0; i < qnu0; ++i) {
				s += Math.abs(q[i] - usrcat[o + i]);
				if (s >= scores[k]) break;
			}
			if (s < scores[k]) {
				scores[k] = s;
				cnfids[k] = j;
			}
		}
	}

	// Sort the scores of compounds of the current chunk.
	const scase = new Uint32Array(new SharedArrayBuffer(Uint32Array.BYTES_PER_ELEMENT * (chunk_end - chunk_beg)));
	for (let k = 0; k < scase.length; ++k) scase[k] = chunk_beg + k;
	scase.sort((val0, val1) => {
		return scores[val0] - scores[val1];
	});

	// Copy the indexes of top hits of the current chunk to a global vector for final sorting.
	zcase.set(scase.slice(0, Math.min(num_hits, chunk_end - chunk_beg)), num_hits * l);
}
