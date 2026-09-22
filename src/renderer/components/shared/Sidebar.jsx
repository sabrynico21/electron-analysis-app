import { NavLink } from 'react-router-dom'
import styles from './Sidebar.module.css'

const NAV_ITEMS = [
  { to: '/',         label: 'Create new analysis', icon: '🔬' },
  { to: '/history',  label: 'History',             icon: '📋' },
  { to: '/settings', label: 'Settings',            icon: '⚙️' },
]

export default function Sidebar() {
  return (
    <nav className={styles.sidebar}>
      <div className={styles.logo}>
        <span>AnalysisApp</span>
      </div>
      <ul className={styles.nav}>
        {NAV_ITEMS.map((item) => (
          <li key={item.to}>
            <NavLink
              to={item.to}
              end={item.to === '/'}
              className={({ isActive }) => `${styles.link} ${isActive ? styles.active : ''}`}
            >
              <span className={styles.icon}>{item.icon}</span>
              <span>{item.label}</span>
            </NavLink>
          </li>
        ))}
      </ul>
    </nav>
  )
}
