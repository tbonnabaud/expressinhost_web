import { createRouter, createWebHistory } from 'vue-router'
import HomeView from '@/views/HomeView.vue'
import CodonTableView from '@/views/CodonTableView.vue'
import TuningView from '@/views/TuningView.vue'
import ResultListView from '@/views/ResultListView.vue'
import ResultView from '@/views/ResultView.vue'
import RegisterView from '@/views/RegisterView.vue'
import AdminView from '@/views/AdminView.vue'
import ForgottenPasswordView from '@/views/ForgottenPasswordView.vue'
import ResetPasswordView from '@/views/ResetPasswordView.vue'

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  // history: createWebHashHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: '/',
      name: 'home',
      component: HomeView,
    },
    {
      path: '/register',
      name: 'register',
      component: RegisterView,
    },
    {
      path: '/forgotten-password',
      name: 'forgotten-password',
      component: ForgottenPasswordView,
    },
    {
      path: '/reset-password/:token',
      name: 'reset-password',
      component: ResetPasswordView,
      props: true,
    },
    {
      path: '/admin',
      name: 'admin',
      component: AdminView,
    },
    {
      path: '/codon-tables',
      name: 'codon-tables',
      component: CodonTableView,
    },
    {
      path: '/tuning',
      name: 'tuning',
      component: TuningView,
    },
    {
      path: '/results',
      name: 'result-list',
      component: ResultListView,
    },
    {
      path: '/results/:id',
      name: 'result',
      component: ResultView,
      props: true,
    },
    {
      path: '/about',
      name: 'about',
      // route level code-splitting
      // this generates a separate chunk (About.[hash].js) for this route
      // which is lazy-loaded when the route is visited.
      component: () => import('../views/AboutView.vue'),
    },
    {
      path: '/:catchAll(.*)',
      name: 'not-found',
      component: () => import('../views/NotFoundView.vue'),
    },
  ],
})

export default router
