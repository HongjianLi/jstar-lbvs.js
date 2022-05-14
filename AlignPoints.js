class Transform3D {
	/**
	 * Initialize a 4x4 matrix to identity.
	 */
	constructor() {
		this.data = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
	}

	/**
	 * @param {Array<number>} quaternion
	 */
	setRotationFromQuaternion(quaternion) {
		const q00 = quaternion[0] * quaternion[0];
		const q11 = quaternion[1] * quaternion[1];
		const q22 = quaternion[2] * quaternion[2];
		const q33 = quaternion[3] * quaternion[3];
		const sumSq = q00 + q11 + q22 + q33;
		const q01 = 2 * quaternion[0] * quaternion[1];
		const q02 = 2 * quaternion[0] * quaternion[2];
		const q03 = 2 * quaternion[0] * quaternion[3];
		const q12 = 2 * quaternion[1] * quaternion[2];
		const q13 = 2 * quaternion[1] * quaternion[3];
		const q23 = 2 * quaternion[2] * quaternion[3];
		this.data[0] = (q00 + q11 - q22 - q33) / sumSq;
		this.data[1] = (q12 + q03) / sumSq;
		this.data[2] = (q13 - q02) / sumSq;
		this.data[4] = (q12 - q03) / sumSq;
		this.data[5] = (q00 - q11 + q22 - q33) / sumSq;
		this.data[6] = (q23 + q01) / sumSq;
		this.data[8] = (q13 + q02) / sumSq;
		this.data[9] = (q23 - q01) / sumSq;
		this.data[10] = (q00 - q11 - q22 + q33) / sumSq;
	}

	/**
	 * Transform a point.
	 * @param {Array<number>} pt - The point to be transformed.
	 */
	transformPoint(pt) {
		const e0 =
			this.data[0] * pt[0] +
			this.data[1] * pt[1] +
			this.data[2] * pt[2] +
			this.data[3];
		const e1 =
			this.data[4] * pt[0] +
			this.data[5] * pt[1] +
			this.data[6] * pt[2] +
			this.data[7];
		const e2 =
			this.data[8] * pt[0] +
			this.data[9] * pt[1] +
			this.data[10] * pt[2] +
			this.data[11];
		pt[0] = e0;
		pt[1] = e1;
		pt[2] = e2;
	}

	/**
	 * @param {Array<number>} move
	 */
	setTranslation(move) {
		this.data[3] = move[0];
		this.data[7] = move[1];
		this.data[11] = move[2];
		this.data[15] = 1.0;
	}
}

/**
 * Compute an optimal alignment (minimum sum of squared distance) between two sets of points in 3D.
 * @param {Array<Array<number>>} refPoints - Reference points.
 * @param {Array<Array<number>>} prbPoints - Points to be aligned to the reference points.
 * @returns {Transform3D} - A transformation.
 */
export default function AlignPoints(refPoints, prbPoints) {
	const npt = refPoints.length;
	const covMat = [
		[0, 0, 0],
		[0, 0, 0],
		[0, 0, 0],
	];
	for (let i = 0; i < npt; ++i) {
		const rpt = refPoints[i];
		const ppt = prbPoints[i];
		covMat[0][0] += ppt[0] * rpt[0];
		covMat[0][1] += ppt[0] * rpt[1];
		covMat[0][2] += ppt[0] * rpt[2];
		covMat[1][0] += ppt[1] * rpt[0];
		covMat[1][1] += ppt[1] * rpt[1];
		covMat[1][2] += ppt[1] * rpt[2];
		covMat[2][0] += ppt[2] * rpt[0];
		covMat[2][1] += ppt[2] * rpt[1];
		covMat[2][2] += ppt[2] * rpt[2];
	}
	const [rptSum, pptSum] = [refPoints, prbPoints].map((points) =>
		points.reduce(
			(sum, pt) => {
				sum[0] += pt[0];
				sum[1] += pt[1];
				sum[2] += pt[2];
				return sum;
			},
			[0, 0, 0]
		)
	);
	const quad = [
		[0, 0, 0, 0],
		[0, 0, 0, 0],
		[0, 0, 0, 0],
		[0, 0, 0, 0],
	];
	let temp;
	temp = pptSum[0] / npt;
	const PxRx = covMat[0][0] - temp * rptSum[0];
	const PxRy = covMat[0][1] - temp * rptSum[1];
	const PxRz = covMat[0][2] - temp * rptSum[2];
	temp = pptSum[1] / npt;
	const PyRx = covMat[1][0] - temp * rptSum[0];
	const PyRy = covMat[1][1] - temp * rptSum[1];
	const PyRz = covMat[1][2] - temp * rptSum[2];
	temp = pptSum[2] / npt;
	const PzRx = covMat[2][0] - temp * rptSum[0];
	const PzRy = covMat[2][1] - temp * rptSum[1];
	const PzRz = covMat[2][2] - temp * rptSum[2];
	quad[0][0] = -2.0 * (PxRx + PyRy + PzRz);
	quad[1][1] = -2.0 * (PxRx - PyRy - PzRz);
	quad[2][2] = -2.0 * (PyRy - PzRz - PxRx);
	quad[3][3] = -2.0 * (PzRz - PxRx - PyRy);
	quad[0][1] = quad[1][0] =  2.0 * (PyRz - PzRy);
	quad[0][2] = quad[2][0] =  2.0 * (PzRx - PxRz);
	quad[0][3] = quad[3][0] =  2.0 * (PxRy - PyRx);
	quad[1][2] = quad[2][1] = -2.0 * (PxRy + PyRx);
	quad[1][3] = quad[3][1] = -2.0 * (PzRx + PxRz);
	quad[2][3] = quad[3][2] = -2.0 * (PyRz + PzRy);
	const eigenVals = [quad[0][0], quad[1][1], quad[2][2], quad[3][3]];
	const eigenVecs = [
		[1, 0, 0, 0],
		[0, 1, 0, 0],
		[0, 0, 1, 0],
		[0, 0, 0, 1],
	];
	let offDiagNorm, diagNorm;
	let b, dma, q, t, c, s;
	let atemp, vtemp, dtemp;
	let i, j, k;
	for (let l = 0; l < 50; l++) {
		diagNorm = 0.0;
		offDiagNorm = 0.0;
		for (j = 0; j <= 3; ++j) {
			diagNorm += Math.abs(eigenVals[j]);
			for (i = 0; i <= j - 1; ++i) {
				offDiagNorm += Math.abs(quad[i][j]);
			}
		}
		if (offDiagNorm / diagNorm <= 1e-6) {
			break;
		}
		for (j = 1; j <= 3; ++j) {
			for (i = 0; i <= j - 1; ++i) {
				b = quad[i][j];
				if (Math.abs(b) > 0.0) {
					dma = eigenVals[j] - eigenVals[i];
					if (Math.abs(dma) + Math.abs(b) <= Math.abs(dma)) {
						t = b / dma;
					} else {
						q = (0.5 * dma) / b;
						t = 1.0 / (Math.abs(q) + Math.sqrt(1.0 + q * q));
						if (q < 0.0) {
							t = -t;
						}
					}
					c = 1.0 / Math.sqrt(t * t + 1.0);
					s = t * c;
					quad[i][j] = 0.0;
					for (k = 0; k <= i - 1; ++k) {
						atemp = c * quad[k][i] - s * quad[k][j];
						quad[k][j] = s * quad[k][i] + c * quad[k][j];
						quad[k][i] = atemp;
					}
					for (k = i + 1; k <= j - 1; ++k) {
						atemp = c * quad[i][k] - s * quad[k][j];
						quad[k][j] = s * quad[i][k] + c * quad[k][j];
						quad[i][k] = atemp;
					}
					for (k = j + 1; k <= 3; ++k) {
						atemp = c * quad[i][k] - s * quad[j][k];
						quad[j][k] = s * quad[i][k] + c * quad[j][k];
						quad[i][k] = atemp;
					}
					for (k = 0; k <= 3; ++k) {
						vtemp = c * eigenVecs[k][i] - s * eigenVecs[k][j];
						eigenVecs[k][j] =
							s * eigenVecs[k][i] + c * eigenVecs[k][j];
						eigenVecs[k][i] = vtemp;
					}
					dtemp =
						c * c * eigenVals[i] +
						s * s * eigenVals[j] -
						2.0 * c * s * b;
					eigenVals[j] =
						s * s * eigenVals[i] +
						c * c * eigenVals[j] +
						2.0 * c * s * b;
					eigenVals[i] = dtemp;
				}
			}
		}
	}
	for (j = 0; j <= 2; ++j) {
		k = j;
		dtemp = eigenVals[k];
		for (i = j + 1; i <= 3; ++i) {
			if (eigenVals[i] < dtemp) {
				k = i;
				dtemp = eigenVals[k];
			}
		}
		if (k > j) {
			eigenVals[k] = eigenVals[j];
			eigenVals[j] = dtemp;
			for (i = 0; i <= 3; ++i) {
				dtemp = eigenVecs[i][k];
				eigenVecs[i][k] = eigenVecs[i][j];
				eigenVecs[i][j] = dtemp;
			}
		}
	}
	const trans = new Transform3D();
	trans.setRotationFromQuaternion([
		eigenVecs[0][0],
		eigenVecs[1][0],
		eigenVecs[2][0],
		eigenVecs[3][0],
	]);
	trans.transformPoint(pptSum);
	trans.setTranslation([
		(rptSum[0] - pptSum[0]) / npt,
		(rptSum[1] - pptSum[1]) / npt,
		(rptSum[2] - pptSum[2]) / npt,
	]);
	return trans;
}
