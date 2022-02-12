'use strict'

const N_MAX_SEARCH = 25
/**
* @description сумма векторов
*/
const summV = function(U, V) {
	const nV = U.length
	const result = []
	for(let i = 0; i < nV; i++) {
		result.push(U[i] + V[i])
	}	
	return result
}
/**
* @description абсолютная разность двух векторов
*/ 
const deltaV = function(U, V) {
	const dX = U[0] - V[0], dY = U[1] - V[1], dZ = U[2] - V[2]
	return Math.sqrt(dX * dX + dY * dY + dZ * dZ)
}
/**
* @description умножение вектора на скаляр
*/ 
const multV = function(U, k) {
	const nV = U.length
	const result = []
	for(let i = 0; i < nV; i++) {
		result.push(U[i] * k)
	}	
	return result	
}
/**
* @description вектор, каждый компонент которого - произведение компонентов входных векторов
*/ 
const multV2V = function(U, V) {
	const nV = U.length
	const result = []
	for(let i = 0; i < nV; i++) {
		result.push(U[i] * V[i])
	}	
	return result	
}
/**
* @description обновить значение вектора
*/ 
const updV = function(U, V) {
	const nV = U.length
	for(let i = 0; i < nV; i++) {
		U[i] += V[i]
	}	
	return U	
}
// Получение градиента
const gradient = function(func, V, eps) {
	const f0 = func(V)
	const result = []
	const nVar = V.length
	for(let i = 0; i < nVar; i++) {
		V[i] += eps
		result.push((func(V) - f0) / eps)
		V[i] -= eps
	}
	
	return result
}
/**
* @description Одномерный поиск на отрезке методом дихотомии
* @param {Function} func исследуемая функция (от вектора действительных чисел)
* @param {Number} x0 начало отрезка поиска
* @param {Number} x1 конец отрезка поиска
* @param {Number} eps оценка погрещности при вычислениях
*/
const dychoSolver_1d = function(func, x0, x1, eps) {
	let x_0 = x0
	let x_1 = x1
	let x_05 = 0.5 * (x_0 + x_1)
	while(Math.abs(x_0 - x_1) > eps) {
		const f0 = func(x_05)
		const f1 = func(x_05 + eps)
		if(f1 > f0) {
			x_0 = x_05
		} else {
			x_1 = x_05
		}
		x_05 = 0.5 *(x_0 + x_1)
	}
	return x_05
}
/**
* @description метод градиентного поиска с подбором шага по каждой из координат методом дихотомии
* @param {Function} func исследуемая функция (от вектора действительных чисел)
* @param {Array<Number>} V0 начальный вектор аргументов
* @param {Array<Number>} delta0 минимальный шаг по каждой координате
* @param {Array<Number>} delta1 максимальный шаг по каждой координате
* @param {Number} eps оценка погрещности при вычислениях
* @param {nMaxGrad} максимальное число шагов при градиентном поиске
*/
const gradientStepDycho = function(func, V0, delta0, delta1, epsDycho, epsGrad, nMaxGrad) {
	const result = [{
		step: 0,
		V: V0.slice(),
		val: func(V0)
	}]
	let i = 0
	const V = V0.slice()
	const nCoord = V0.length
	
	while(i++ < nMaxGrad) {
		const grad = gradient(func, V, epsGrad)
		const delta_active = delta0.slice()
		for(let j = 0; j < nCoord; j++) {
			const stepValue = function(R) {
				delta_active[j] = R
				const V_1 = summV(V, multV2V(grad, delta_active))
				return func(V_1)
			}
			delta_active[j] = dychoSolver_1d(stepValue, delta0[j], delta1[j], epsDycho)
		}

		updV(V, multV2V(grad, delta_active))

		result.push({
			step: i,
			V: V.slice(),
			val: func(V)
		})

		const dF = Math.abs(result[i].val - result[i - 1].val)
		if(dF < epsGrad) { break }
	}
	return result
}


const testFunc = x => -75 + 7.5 * x - 0.0125 * (x**4)
const testFunc_XY = V => {
	const X = V[0], Y = V[1]
	return -75 + 12.5 * X - 0.05 * X * X + 3.5 * Y - 0.0125 * Y * Y - 0.15 * X * Y
}

const sqrtApprox = (koef, X) => koef[0] + koef[1] * X + koef[2] * X * X

const regressError = V => koef => {
	const n = V.length
	let delta = 0
	for(let i = 0; i < n; i++) {
		delta += Math.abs(V[i][1] - sqrtApprox(koef, V[i][0]))
	}
	return -delta
}

const testSet = [
	[-15, 2],
	[-13, -7],
	[-11, -9],
	[-9, -6],
	[-7, -5],
	[-5, -7],
	[-3, -3],
	[-1, -3],
	[1, -2],
	[3, 4],
	[5, 1],
	[7, 5],
	[8, 4],
	[11, 7],
]

const regressTest = regressError(testSet)
const test_opt = gradientStepDycho(
	regressTest,
	[5, 5, 5],
	[0.0, 0.0, 0.0],
	[10, 10, 10],
	1E-5,
	1E-3,
	N_MAX_SEARCH
)
console.log(test_opt)

const regressStat = (V, X, Y) => {
	const nV = V.length
	for(let i = 0; i < nV; i++) {
		
	}
}