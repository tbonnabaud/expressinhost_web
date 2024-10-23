import axios, { type AxiosRequestConfig } from 'axios'
import type { RunTrainingForm, TuningOutput } from './interfaces'

const client = axios.create({
  baseURL: '/api',
})

/**
 * The function `request` sends an asynchronous HTTP request using Axios and handles errors
 * appropriately.
 * @param {AxiosRequestConfig} config - The `config` parameter in the `request` function is of type
 * `AxiosRequestConfig`, which is an interface defining the configuration for an Axios HTTP request. It
 * typically includes properties like `url`, `method`, `headers`, `params`, `data`, etc., that are used
 * to customize
 * @returns The `request` function returns the data from the response if the request is successful. If
 * there is an error, it checks if the error is an Axios error with a response, and if so, it logs the
 * error message and response data. If it's a network error or the request failed, it logs a generic
 * error message. In both cases, it returns `null`.
 */
async function request(config: AxiosRequestConfig) {
  try {
    const response = await client.request(config)
    return response.data
  } catch (error) {
    if (axios.isAxiosError(error) && error.response) {
      console.error(error.message, error.response.data)
    } else {
      console.error('Network error or request failed')
    }
    return null
  }
}

const REQUESTS = {
  get: async (url: string, params?: object) =>
    request({ method: 'get', url, params }),
  post: async (url: string, data?: object) =>
    request({ method: 'post', url, data }),
  put: async (url: string, data?: object) =>
    request({ method: 'put', url, data }),
  delete: async (url: string) => request({ method: 'delete', url }),
}

export const API = {
  runTraining: async (form: RunTrainingForm) =>
    (await REQUESTS.post('/run-tuning', form)) as TuningOutput,
}
