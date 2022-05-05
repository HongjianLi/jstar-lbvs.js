export default function local_time_string() {
	const date = new Date();
	return `${[ date.getFullYear(), date.getMonth() + 1, date.getDate() ].map((component) => {
		return component.toString().padStart(2, '0');
	}).join('-')} ${[ date.getHours(), date.getMinutes(), date.getSeconds() ].map((component) => {
		return component.toString().padStart(2, '0');
	}).join(':')}.${date.getMilliseconds().toString().padStart(3, '0')}`;
}
