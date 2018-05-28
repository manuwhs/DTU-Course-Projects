export class HelloWorld extends React.Component {
  render() {
    return (
      <div className="helloworld">
        Hello {this.props.name}
      </div>
    );
  }
}