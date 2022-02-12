'use strict'
/**
* @description максимальное число итераций градиентного метода
*/
const N_MAX_SEARCH = 75
/**
* @description сумма векторов
* @param U {Array<Number>}
* @param V {Array<Number>}
* @returns {Array<Number>}
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
* @param U {Array<Number>}
* @param V {Array<Number>}
* @returns {Array<Number>}
*/ 
const deltaV = function(U, V) {
	const dX = U[0] - V[0], dY = U[1] - V[1], dZ = U[2] - V[2]
	return Math.sqrt(dX * dX + dY * dY + dZ * dZ)
}
/**
* @description умножение вектора на скаляр
* @param U {Array<Number>}
* @param k {Number}
* @returns {Array<Number>}
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
* @description вектор, каждый компонент которого - произведение компонентов входных векторов U,V -> [U[i]*V[i], i = 1..N]
* @param U {Array<Number>}
* @param V {Array<Number>}
* @returns {Array<Number>}
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
* @param U {Array<Number>}
* @param V {Array<Number>}
* @returns {Array<Number>}
*/ 
const updV = function(U, V) {
	const nV = U.length
	for(let i = 0; i < nV; i++) {
		U[i] += V[i]
	}	
	return U	
}
/**
* @description Получить градиент функции от вектора действительных чисел
* @param func {Function(Array<Number>)} функция, для которой получен градиент
* @param V {Array<Number>} точка, в которой вычисляется градиент
* @param eps {Number} шаг разности
* @returns {Array<Number>}
*/
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
* @param {Function(Array<Nymber>)} func исследуемая функция (от вектора действительных чисел)
* @param {Number} x0 начало отрезка поиска
* @param {Number} x1 конец отрезка поиска
* @param {Number} eps оценка погрещности при вычислениях
* @return {Number}
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
* @param {Function(Array<Number>)} func исследуемая функция (от вектора действительных чисел)
* @param {Array<Number>} V0 начальный вектор аргументов
* @param {Array<Number>} delta0 минимальный шаг по каждой координате
* @param {Array<Number>} delta1 максимальный шаг по каждой координате
* @param {Number} epsDycho оценка погрещности при выборе шага смещения
* @param {Number} epsDycho оценка погрещности при вычислени градиента
* @param {nMaxGrad} максимальное число шагов при градиентном поиске
* @return {Array<Object.{Number, Number, Array<Number>}>}
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
	let maxVal = -1E10
	
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

		const activeVal = func(V)
		result.push({
			step: i,
			V: V.slice(),
			val: activeVal
		})
		
		const prevMax = maxVal
		maxVal = Math.max(maxVal, activeVal)

		if(Math.abs(maxVal - prevMax) < epsGrad) { break }
	}
	return result
}

/**
* @description квадратическая функция
*/
const sqrtApprox = (koef, X) => koef[0] + koef[1] * X + koef[2] * X * X

/**
* @descripotion ошибка квадратической регрессии как функция от начального набора данных и коэффициентов регрессии
*/
const regressError = V => koef => {
	const n = V.length
	let delta = 0
	for(let i = 0; i < n; i++) {
		delta += Math.abs(V[i][1] - sqrtApprox(koef, V[i][0]))
	}
	return -delta
}
/**
* @description тестовый набор данных
*/
const testSet = [
	[-15, 2],
	[-13, -7],
	[-11, -9],
	[-9, -6],
	[-7, -5],
	[-5, -7],
	[-4, -7],
	[-3, -3],
	[-1, -3],
	[1, -2],
	[3, 4],
	[5, 1],
	[7, 5],
	[8, 4],
	[11, 7],
	[12, 6],
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