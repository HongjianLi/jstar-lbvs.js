import os from 'os';
import path from 'path';
import { promises as fs } from 'fs';
import { program } from 'commander';
import { MongoClient } from "mongodb";
import initRDKitModule from '@rdkit/rdkit';
import { Piscina } from 'piscina';
import { CompoundDatabase } from './CompoundDatabase.js';
import local_time_string from './utility.js';
import AlignPoints from './AlignPoints.js';

/**
 * Calculate the square distance between two points.
 * @param {Array<number>} p0 
 * @param {Array<number>} p1 
 * @returns {number} Square distance.
 */
function dist2(p0, p1) {
	const d0 = p0[0] - p1[0];
	const d1 = p0[1] - p1[1];
	const d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

/**
 * Calculate four reference points of a molecule.
 * @param {Array<Array<number>>} cnf - Conformer.
 * @param {Array<number>} heavyAtoms - Indexes of heavy atoms.
 * @returns {Array<Array<number>>} - Four reference points.
 */
function calcRefPoints(conf, heavyAtoms) {
	const num_points = heavyAtoms.length;
	let ctd = Array(3).fill(0);
	let cst = Array(3).fill(0);
	let fct = Array(3).fill(0);
	let ftf = Array(3).fill(0);
	for (const i of heavyAtoms) {
		const a = conf[i];
		ctd[0] += a[0];
		ctd[1] += a[1];
		ctd[2] += a[2];
	}
	ctd[0] /= num_points;
	ctd[1] /= num_points;
	ctd[2] /= num_points;
	let cst_dist = Infinity;
	let fct_dist = -Infinity;
	let ftf_dist = -Infinity;
	for (const i of heavyAtoms) {
		const a = conf[i];
		const this_dist = dist2(a, ctd);
		if (this_dist < cst_dist) {
			cst = a;
			cst_dist = this_dist;
		}
		if (this_dist > fct_dist) {
			fct = a;
			fct_dist = this_dist;
		}
	}
	for (const i of heavyAtoms) {
		const a = conf[i];
		const this_dist = dist2(a, fct);
		if (this_dist > ftf_dist) {
			ftf = a;
			ftf_dist = this_dist;
		}
	}
	return [ ctd, cst, fct, ftf ];
}

// Process program options.
const options = program
	.name('jstar-lbvs.js')
	.description("jstar's daemon for ligand-based virtual screening (LBVS), written in JavaScript.")
	.version('1.0.0')
	.option('--databases <string>', 'Path to the compound database directory', '../jstar/databases')
	.option('--host <string>', 'DBMS host', 'localhost')
	.option('--port <number>', 'DBMS port', cur => parseInt(cur), 27017)
	.option('--user <string>', 'DBMS user', 'jstard')
	.option('--pass <string>', 'DBMS password')
	.option('--threads <number>', 'number of worker threads to use', cur => parseInt(cur), os.cpus().length)
	.parse()
	.opts();

// Initialize constants.
console.log(`${local_time_string()} Initializing`);
const num_usrs = 2;
const usr_names = [ "USR", "USRCAT" ];
const qn = [ 12, 60 ];
const qv = qn.map(n => 1.0 / n);
const num_refPoints = 4;
const SubsetSMARTS = [
	"[!#1]", // heavy
	"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
	"[a]", // aromatic
	"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
	"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
];
const num_subsets = SubsetSMARTS.length;
const num_hits = 100;

// Wrap SMARTS strings to Mol objects.
const rdkit = await initRDKitModule();
const SubsetMols = SubsetSMARTS.map(smarts => rdkit.get_qmol(smarts)); // https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/common.h#L173

// Read compound database directory.
const databases = await fs.readdir(options.databases).then(subDirs => subDirs.map(subDir => new CompoundDatabase(path.join(options.databases, subDir))));
for (let k = 0; k < databases.length; ++k) { await databases[k].read_descriptors(); }

// Connect to mongodb and authenticate user.
console.log(`${local_time_string()} Connecting to ${options.host}:${options.port} and authenticating ${options.user}`);
const mongoClient = new MongoClient(`mongodb://${options.host}:${options.port}/?authSource=jstar&maxPoolSize=3`); // https://www.mongodb.com/docs/drivers/node/current/fundamentals/connection/ Always URI encode the username and password using the encodeURIComponent method to ensure they are correctly parsed.
await mongoClient.connect();
const jstar = mongoClient.db('jstar');
const coll = jstar.collection('lbvs');

const jobid_filter = { 
	startDate : {
		$exists : false,
	},
	database: {
		$in: databases.map(db => db.name),
	}
};
const jobid_foau_options = { // By default, the original document is returned
	$sort: {
		submitDate: 1,
	},
	$projection: {
		_id: 1,
		qryMolSdf: 1,
		database: 1,
		score: 1,
	},
};

// Initialize variables.
/**
 * array<vector<int>, num_subsets> subsets;
 */
const subsets = Array(num_subsets);
/**
 * array<vector<double>, num_refPoints> dista;
 */
const dista = Array(num_refPoints);
/**
 * alignas(32) array<double, 60> q;
 */
const q = new Float32Array(new SharedArrayBuffer(Float32Array.BYTES_PER_ELEMENT * 60)); // The usrcat feature vector of the query molecule.

// Create worker threads for later use.
const piscina = new Piscina({
	filename: new URL('worker.js', import.meta.url).href,
	maxThreads: options.threads,
});

// Enter event loop.
console.log(`${local_time_string()} Entering event loop`);
let sleeping = false;
while (true) {
	// Fetch an incompleted job in a first-come-first-served manner.
	if (!sleeping) console.log(`${local_time_string()} Fetching an incompleted job`);
	const startDate = new Date();
	const jobid_document = await coll.findOneAndUpdate(jobid_filter, { $set: { startDate } }, jobid_foau_options); // https://mongodb.github.io/node-mongodb-native/4.5/classes/Collection.html#findOneAndUpdate
	if (!jobid_document.value) {
		// No incompleted jobs. Sleep for a while.
		if (!sleeping) console.log(`${local_time_string()} Sleeping`);
		sleeping = true;
		await new Promise(resolve => setTimeout(resolve, 2000));
		continue;
	}
	sleeping = false;
	const jobid_view = jobid_document.value;

	// Obtain job properties.
	const _id = jobid_view["_id"];
	console.log(`${local_time_string()} Executing job ${_id}`);
	const qryMolSdf = jobid_view["qryMolSdf"];
	const cpdb_name = jobid_view["database"];
	const score = jobid_view["score"];
	console.assert(usr_names.includes(score));
	const usr0 = usr_names.indexOf(score); // Specify the primary sorting score. 0: USR; 1: USRCAT.
	const usr1 = usr0 ^ 1;
	const qnu0 = qn[usr0];
	const qnu1 = qn[usr1];

	// Obtain a constant reference to the selected database.
	console.log(`${local_time_string()} Finding the selected compound database`);
	const cpdb = databases.find(cpdb => cpdb.name === cpdb_name);
	const usrcat = cpdb['usrcat.f32'];

	// Read the user-supplied SDF file.
	console.log(`${local_time_string()} Reading the query file`);
	const qrySdfArr = qryMolSdf.split(/\$\$\$\$\r?\n/).filter(s => s.length).map(s => s + '$$$$\n');

	// Initialize vectors to store compounds' primary score and their corresponding conformer.
	/**
	 * Primary score of compounds.
	 * vector<double> scores(cpdb.num_compounds);
	 * @type {Array<number>}
	 */
	const scores = new Float32Array(new SharedArrayBuffer(Float32Array.BYTES_PER_ELEMENT * cpdb.num_compounds));
	/**
	 * ID of conformer with the best primary score.
	 * vector<size_t> cnfids(cpdb.num_compounds);
	 * @type {Array<number>}
	 */
	const cnfids = new Uint32Array(new SharedArrayBuffer(Uint32Array.BYTES_PER_ELEMENT * cpdb.num_compounds));

	// Initialize the number of chunks and the number of compounds per chunk.
	const chunk_size = Math.max(1 + Math.floor((cpdb.num_compounds - 1) / (options.threads << 2)), num_hits);
	const num_chunks = 1 + Math.floor((cpdb.num_compounds - 1) / chunk_size);
	console.assert(chunk_size * num_chunks >= cpdb.num_compounds);
	console.assert(chunk_size >= num_hits);
	console.log(`${local_time_string()} Using ${num_chunks} chunks and a chunk size of ${chunk_size}`);
	/**
	 * The last chunk might have fewer than num_hits records.
	 * vector<size_t> zcase(num_hits * (num_chunks - 1) + Math.min(num_hits, cpdb.num_compounds - chunk_size * (num_chunks - 1)));
	 */
	const zcase = new Uint32Array(new SharedArrayBuffer(Uint32Array.BYTES_PER_ELEMENT * (num_hits * (num_chunks - 1) + Math.min(num_hits, cpdb.num_compounds - chunk_size * (num_chunks - 1))))); // The last chunk might have fewer than num_hits records.

	// Process each of the query compounds sequentially.
	let hitMolSdf = '';
	const num_qry_mols = qrySdfArr.length; // num_qry_mols is the number of query molecules submitted by the user. These query molecules are not necessarily all processed, given the limitation of maximum 16MB MongoDB document size of the result.
	let query_number = 0;
	while (query_number < num_qry_mols) {
		console.log(`${local_time_string()} Parsing query compound ${query_number}`);
		const qryMol = rdkit.get_mol(qrySdfArr[query_number]); // get_mol has a different implementation from C++ SDMolSupplier(sanitize=true, removeHs=true, strictParsing=true). So the Mol object returned by get_mol is not totally identical to the ROMol instance returned by SDMolSupplier.next(), thus they do not return identical Morgan fingerprints. https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/common.h#L89
		const qryCnf = JSON.parse(qryMol.get_json()).molecules[0].conformers[0].coords;
		const qryDes = JSON.parse(qryMol.get_descriptors());

		// Get the number of atoms, including and excluding hydrogens.
		const num_atoms = qryDes.NumAtoms;
		const num_heavy_atoms = qryDes.NumHeavyAtoms;
		console.assert(num_atoms === qryCnf.length);
		console.assert(num_heavy_atoms);
		console.log(`${local_time_string()} Found ${num_atoms} atoms and ${num_heavy_atoms} heavy atoms`);

		// Calculate Morgan fingerprint.
		console.log(`${local_time_string()} Calculating Morgan fingerprint`);
		const qryFp = qryMol.get_morgan_fp(2, 2048); // https://github.com/rdkit/rdkit/blob/master/Code/MinimalLib/minilib.cpp#L148

		// Classify atoms to pharmacophoric subsets.
		console.log(`${local_time_string()} Classifying atoms into subsets`);
		for (let k = 0; k < num_subsets; ++k) {
			const matchVect = JSON.parse(qryMol.get_substruct_matches(SubsetMols[k])); // get_substruct_matches returns a string, e.g. [{"atoms":[0],"bonds":[]},{"atoms":[1],"bonds":[]}]
			subsets[k] = matchVect.map(m => m.atoms[0]);
			console.log(`${local_time_string()} Found ${matchVect.length} atoms for subset ${k}`);
		}
		const subset0 = subsets[0];
		console.assert(subset0.length == num_heavy_atoms);

		// Calculate the four reference points.
		console.log(`${local_time_string()} Calculating ${num_refPoints} reference points`);
		const qryRefPoints = calcRefPoints(qryCnf, subset0);

		// Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
		console.log(`${local_time_string()} Calculating ${num_heavy_atoms * num_refPoints} pairwise distances`);
		for (let k = 0; k < num_refPoints; ++k) {
			const refPoint = qryRefPoints[k];
			const distp = dista[k] = Array(num_atoms);
			for (let i = 0; i < num_heavy_atoms; ++i) {
				distp[subset0[i]] = Math.sqrt(dist2(qryCnf[subset0[i]], refPoint));
			}
		}

		// Loop over pharmacophoric subsets and reference points.
		console.log(`${local_time_string()} Calculating ${3 * num_refPoints * num_subsets} moments of USRCAT feature`);
		let qo = 0;
		for (const subset of subsets) {
			const n = subset.length;
			for (let k = 0; k < num_refPoints; ++k) {
				// Load distances from precalculated ones.
				const distp = dista[k];
				const dists = Array(n);
				for (let i = 0; i < n; ++i) {
					dists[i] = distp[subset[i]];
				}

				// Compute moments.
				const m = Array(3).fill(0);
				if (n > 2) {
					const v = 1.0 / n;
					for (let i = 0; i < n; ++i) {
						const d = dists[i];
						m[0] += d;
					}
					m[0] *= v;
					for (let i = 0; i < n; ++i) {
						const d = dists[i] - m[0];
						m[1] += d * d;
					}
					m[1] = Math.sqrt(m[1] * v);
					for (let i = 0; i < n; ++i) {
						const d = dists[i] - m[0];
						m[2] += d * d * d;
					}
					m[2] = Math.cbrt(m[2] * v);
				} else if (n == 2) {
					m[0] = 0.5 * (dists[0] + dists[1]);
					m[1] = 0.5 * Math.abs(dists[0] - dists[1]);
				} else if (n == 1) {
					m[0] = dists[0];
				}
				for (const e of m) {
					q[qo++] = e;
				}
			}
		}
		console.assert(qo == qn[1]);

		// Compute USR and USRCAT scores.
		console.log(`${local_time_string()} Screening ${cpdb.name} and calculating ${cpdb.num_compounds} ${usr_names[usr0]} scores from ${cpdb.num_conformers} conformers`);
		scores.fill(Infinity);
		await Promise.all([...Array(num_chunks).keys()].map(l => {
			return piscina.run({ chunk_size, l, num_compounds: cpdb.num_compounds, scores, usrcat, qnu0, q, cnfids, zcase, num_hits }, { name: 'calculate' });
		}));

		// Sort the top hits from chunks.
		console.log(`${local_time_string()} Sorting ${zcase.length} hits by ${usr_names[usr0]} score`);
		zcase.sort((val0, val1) => {
			return scores[val0] - scores[val1];
		});

		// Write hit molecules to a string for output.
		console.log(`${local_time_string()} Writing hit molecules to a string`);
		let hitMolSdfPerQry = '';
		for (let l = 0; l < num_hits; ++l) {
			// Obtain indexes to the hit compound and the hit conformer.
			const k = zcase[l];
			const j = cnfids[k];

			// Calculate the secondary score of the saved conformer, which has the best primary score.
			const o = 60 * j;
			let s = 0;
			for (let i = 0; i < qnu1; ++i) {
				s += Math.abs(q[i] - usrcat[o + i]);
			}
			const u0score = 1 / (1 + scores[k] * qv[usr0]); // Primary score of the current compound.
			const u1score = 1 / (1 + s * qv[usr1]); // Secondary score of the current compound.

			// Read SDF content of the hit conformer.
			const hitSdf = await cpdb.read_conformer(j);

			// Construct a RDKit ROMol object.
			const hitMol = rdkit.get_mol(hitSdf);
			const hitCnf = JSON.parse(hitMol.get_json()).molecules[0].conformers[0].coords;
			const hitDes = JSON.parse(hitMol.get_descriptors());

			// Calculate Morgan fingerprint.
			const hitFp = hitMol.get_morgan_fp(2, 2048);
			console.assert(hitFp.length === qryFp.length);

			// Calculate Tanimoto similarity.
			let qryFpNumOnBits = 0, hitFpNumOnBits = 0, bothFpNumOnBits = 0; // https://www.rdkit.org/docs/cppapi/BitOps_8h.html#aee73e04c05be3284829cfe7b16beb72f
			for (let i = 0; i < qryFp.length; ++i) {
				const qryFpOnBit = qryFp[i] === '1';
				const hitFpOnBit = hitFp[i] === '1';
				qryFpNumOnBits += (qryFpOnBit);
				hitFpNumOnBits += (hitFpOnBit);
				bothFpNumOnBits += (qryFpOnBit && hitFpOnBit);
			}
			const ts = bothFpNumOnBits / (qryFpNumOnBits + hitFpNumOnBits - bothFpNumOnBits);

			// Remove hydrogens to calculate canonical SMILES and descriptors.
			const hitMolNoH = rdkit.get_mol(hitMol.remove_hs());

/*			// Calculate canonical SMILES, molecular formula and descriptors. This can be done either by calculating them on the fly using the molecule with hydrogens removed, or by reading the precalculated values from *.u16 and *.f32 files.
			hitMol.setProp<unsigned int>("query", query_number);
			hitMol.setProp<double>("usrScore", usr1 ? u0score : u1score);
			hitMol.setProp<double>("usrcatScore", usr1 ? u1score : u0score);
			hitMol.setProp<double>("tanimotoScore", ts);
			hitMol.setProp<string>("database", cpdb.name);
			hitMol.setProp<string>("canonicalSMILES", hitMolNoH.get_smiles()); // Default parameters are: const ROMol& mol, bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1, bool canonical = true, bool allBondsExplicit = false, bool allHsExplicit = false, bool doRandom = false. https://www.rdkit.org/docs/cppapi/namespaceRDKit.html#a3636828cca83a233d7816f3652a9eb6b
			hitMol.setProp<string>("molFormula", calcMolFormula(hitMol)); // Calculate hydrogens in the molecular formula.
			hitMol.setProp<unsigned int>("numAtoms", cpdb['natm.u16'][k]); // hitDes.NumHeavyAtoms
			hitMol.setProp<unsigned int>("numHBD", cpdb['nhbd.u16'][k]); // hitDes.NumHBD
			hitMol.setProp<unsigned int>("numHBA", cpdb['nhba.u16'][k]); // hitDes.NumHBA
			hitMol.setProp<unsigned int>("numRotatableBonds", cpdb['nrtb.u16'][k]); // cpdb['nrtb.u16'][k] was precalculated from SMILES before adding hydrogens. Adding hydrogens may lead to more rotatable bonds. As a result, cpdb['nrtb.u16'][k] == calcNumRotatableBonds(hitMolNoH) <= calcNumRotatableBonds(hitMol) // hitDes.NumRotatableBonds
			hitMol.setProp<unsigned int>("numRings", cpdb['nrng.u16'][k]); // hitDes.NumRings
			hitMol.setProp<double>("exactMW", cpdb['xmwt.f32'][k]); // hitDes.exactmw
			hitMol.setProp<double>("tPSA", cpdb['tpsa.f32'][k]); // hitDes.tpsa
			hitMol.setProp<double>("clogP", cpdb['clgp.f32'][k]); // hitDes.CrippenClogP */

			// Find heavy atoms.
			const matchVect = JSON.parse(hitMol.get_substruct_matches(SubsetMols[0]));
			const hitHeavyAtoms = matchVect.map(m => m.atoms[0]);
			console.assert(hitHeavyAtoms.length === hitDes.NumHeavyAtoms);
//			for (let i = 0; i < num_matches; ++i) {
//				console.assert(hitHeavyAtoms[i] == i); // Comment this assertion to avoid warning: comparison of integer expressions of different signedness. hitHeavyAtoms can be constructed using iota(hitHeavyAtoms.begin(), hitHeavyAtoms.end(), 0); because for RDKit-generated SDF compounds, heavy atom are always the first few atoms.
//			}

			// Calculate the four reference points.
			const hitRefPoints = calcRefPoints(hitCnf, hitHeavyAtoms);

			// Calculate a 3D transform from the four reference points of the hit conformer to those of the query compound.
			// Ported from rdkit-2020.03.6/Code/Numerics/Alignment/AlignPoints.cpp
			const trans = AlignPoints(qryRefPoints, hitRefPoints);

			// Apply the 3D transform to all atoms of the hit conformer.
			// Ported from rdkit-2020.03.6/Code/GraphMol/MolTransforms/MolTransforms.cpp
			hitCnf.forEach(pt => {
				trans.transformPoint(pt);
			});

			// Write the aligned hit conformer.
			hitMolSdfPerQry += (hitMol.get_molblock() + '$$$$\n').split('\n').map((line, lineIndex) => { // .get_molblock() does not return the trailing $$$$ line.
				const atomIndex = lineIndex - 4;
				if (atomIndex < 0 || atomIndex >= hitCnf.length) return line;
				const point = hitCnf[atomIndex];
				return point.map(e => e.toFixed(4).padStart(10, ' ')).join('').concat(line.slice(30));
			}).join('\n');
		}
		console.log(`${local_time_string()} Wrote ${hitMolSdfPerQry.length} bytes of hit molecules to a string`);

		// If the size of the hitMolSdf field will not exceed 15MB after appending, then allow the appending. Reserve 1MB for the other fields, e.g. qryMolSdf.
		if (hitMolSdf.length + hitMolSdfPerQry.length < 15000000) {
			hitMolSdf += hitMolSdfPerQry;
			console.log(`${local_time_string()} Accumulated ${hitMolSdf.length} bytes of hit molecules for the first ${++query_number} query molecules`);
		} else {
			console.log(`${local_time_string()} Unable to accumulate ${hitMolSdfPerQry.length} bytes to the existing ${hitMolSdf.length} bytes`);
			break;
		}
	}
	const num_qry_mols_processed = query_number; // num_qry_mols_processed is the number of query molecules processed by the daemon.
	console.assert(num_qry_mols_processed <= num_qry_mols);
	console.log(`${local_time_string()} Processed ${num_qry_mols_processed} out of ${num_qry_mols} query molecules`);

	// Update job status.
	console.log(`${local_time_string()} Writing ${hitMolSdf.length} bytes of hit molecules and setting end date`);
	const endDate = new Date();
	const compt_update = await coll.updateOne({ _id }, {
		$set: {
			hitMolSdf,
			endDate,
			numQryMol: num_qry_mols_processed,
			numLibMol: cpdb.num_compounds,
			numLibCnf: cpdb.num_conformers,
		},
	}, {});
	console.assert(compt_update); // {acknowledged: true, modifiedCount: 1, upsertedId: null, upsertedCount: 0, matchedCount: 1}
	console.assert(compt_update.matchedCount === 1);
	console.assert(compt_update.modifiedCount === 1);

	// Calculate runtime in seconds and screening speed in conformers per second.
	const runtime = (endDate.getTime() - startDate.getTime()) * 1e-3; // Convert milliseconds to seconds.
	const speed = cpdb.num_conformers * num_qry_mols_processed / runtime;
	console.log(`${local_time_string()} Completed ${num_qry_mols_processed} ${num_qry_mols_processed == 1 ? "query" : "queries"} in ${runtime.toFixed(3)} seconds`);
	console.log(`${local_time_string()} Screening speed was ${speed.toFixed(0)} conformers per second`);
}
await mongoClient.close();
